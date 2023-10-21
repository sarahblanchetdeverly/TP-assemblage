#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
from networkx import DiGraph, all_simple_paths, lowest_common_ancestor, has_path, random_layout, draw, spring_layout
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Sarah BLANCHET--DEVERLY"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Sarah BLANCHET--DEVERLY"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Sarah BLANCHET--DEVERLY"
__email__ = "sarahblanchetdeverly@icloud.com"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open(fastq_file, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            # Extract the sequence from fastq format
            sequence = lines[i + 1].strip()
            yield sequence
    

def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(len(read) - kmer_size + 1):
        kmer = read[i:i + kmer_size]
        yield kmer


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}
    
    # Iterate through the reads in the FASTQ file
    for read in read_fastq(fastq_file):
        # Iterate through the k-mers in the read
        for kmer in cut_kmer(read, kmer_size):
            # Count the occurrences of each k-mer
            kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
    
    return kmer_dict
    

def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
        # Create a directed graph
    graph = nx.DiGraph()

    # Iterate through the k-mers and their occurrences
    for kmer, occurrence in kmer_dict.items():
        # Split the k-mer into prefix and suffix
        prefix = kmer[:-1]
        suffix = kmer[1:]

        # Check if the nodes for the prefix and suffix exist in the graph
        if not graph.has_node(prefix):
            graph.add_node(prefix)
        if not graph.has_node(suffix):
            graph.add_node(suffix)

        # Add an edge with the weight (occurrence) from prefix to suffix
        graph.add_edge(prefix, suffix, weight=occurrence)

    return graph


def remove_paths(my_graph, path_list, delete_entry_node, delete_sink_node):
    '''
    Deletes all the paths from path_list in my_graph.
    If delete_entry_node is set to True, the first node of each path is deleted.
    If delete_sink_node is set to True, the last node of each path is deleted.
    '''

    for path in path_list:

        edges_list=[(path[i],path[i+1]) for i in range(len(path)-1)]
        my_graph.remove_edges_from(edges_list)

        nodes_list=[path[i] for i in range(1,len(path)-1)]
        my_graph.remove_nodes_from(nodes_list)

        if delete_entry_node:
            my_graph.remove_node(path[0])

        if delete_sink_node:
            my_graph.remove_node(path[-1])

    return my_graph


def select_best_path(my_graph, path_list, path_length_list, path_average_weight_list,
                     delete_entry_node=False, delete_sink_node=False):

    '''
    Selects the best path from a graph list.
    Priority is given to the highest average weight, and then to the longest path length.
    If several paths are equal, a random path is selected.
    '''
    # comparaison of the average weights :
    weight_std = statistics.stdev(path_average_weight_list)

    if weight_std > 0 :

        max_weight = max(path_average_weight_list)
        removed_path_index = []

        for i in range(len(path_average_weight_list)):

            if path_average_weight_list[i] < max_weight:

                my_graph = remove_paths(my_graph, [path_list[i]], delete_entry_node, delete_sink_node)
                removed_path_index.append(i)

        updated_path_list = path_list.copy()
        updated_length_list = path_length_list.copy()
        updated_average_weight_list = path_average_weight_list.copy()

        for index in removed_path_index:

            updated_path_list.remove(path_list[index])
            updated_length_list.remove(path_length_list[index])
            updated_average_weight_list.remove(path_average_weight_list[index])

        path_list = updated_path_list.copy()
        path_length_list = updated_length_list.copy()
        path_average_weight_list = updated_average_weight_list.copy()

    if len(path_list) > 1:

        # comparaison of the lengths :
        length_std = statistics.stdev(path_length_list)

        if length_std > 0:

            max_length = max(path_length_list)
            removed_path_index = []

            for i in range(len(path_list)):

                if path_length_list[i] < max_length:

                    my_graph = remove_paths(my_graph, [path_list[i]],
                                          delete_entry_node,
                                          delete_sink_node)
                    removed_path_index.append(i)

            updated_path_list = path_list.copy()
            updated_length_list = path_length_list.copy()
            updated_average_weight_list = path_average_weight_list.copy()

            for index in removed_path_index:

                updated_path_list.remove(path_list[index])
                updated_length_list.remove(path_length_list[index])
                updated_average_weight_list.remove(path_average_weight_list[index])

            path_list = updated_path_list.copy()
            path_length_list = updated_length_list.copy()
            path_average_weight_list = updated_average_weight_list.copy()

    if len(path_list) > 1: # random selection

        random_index = randint(0, len(path_list) - 1)

        for i in range(len(path_list)):

            if path_list[i] != path_list[random_index]:

                my_graph = remove_paths(my_graph, [path_list[i]], delete_entry_node, delete_sink_node)

        path_list = [path_list[random_index]]
        path_length_list = [path_length_list[random_index]]
        path_average_weight_list = [path_average_weight_list[random_index]]

    return my_graph


def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    '''
    For a given 'bubble' (several paths between two nodes), selects the best path.
    '''
    #finding all the paths between the ancestor and the descendant
    p_list=nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node)
    all_paths=[path for path in p_list]

    if len(all_paths)==1:
        return graph

    len_list=[len(path) for path in all_paths]
    avg_weight_list=[path_average_weight(graph, path) for path in all_paths]

    graph=select_best_path(graph, all_paths, len_list, avg_weight_list)

    return graph


def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """

    nodes = list(graph.nodes())

    number_of_unchecked = 0

    for node in nodes:
        in_edges = [(from_node, node) for from_node in graph.predecessors(node)]
        if len(in_edges) > 1:
            number_of_unchecked += 1

    if number_of_unchecked > 0:
        node_index = 0
        node = nodes[node_index]

        while (len(graph.in_edges(node)) < 2) and (node_index < len(nodes) - 1):
            node_index += 1
            node = nodes[node_index]

        edges = list(graph.in_edges(node))
        combinations = [(edges[i], edges[j]) for i in range(len(edges)) for j in range(i + 1, len(edges))]

        comb = combinations[0]
        node_1 = comb[0][0]
        node_2 = comb[1][0]

        common_ancestor = nx.lowest_common_ancestor(graph, node_1, node_2)
        if common_ancestor is not None:
            graph = solve_bubble(graph, common_ancestor, node)
            number_of_unchecked -= 1

            if number_of_unchecked > 0:
                return simplify_bubbles(graph)
            else:
                return graph
            


def solve_entry_tips(graph,starting_nodes):
    '''
     """Remove entry tips and returns a graph without useless entry paths on the simplify bubble basis

    :param graph: (nx.DiGraph) A directed graph object
    :starting_nodes: (list) A list of nodes without predecessors
    :return: (nx.DiGraph) A directed graph object
    """
    '''

    path_list = [list(nx.all_simple_paths(graph, start, node)) for node in graph.nodes for start in starting_nodes if len(list(graph.predecessors(node))) > 1]
    if path_list:
        path_length = [len(path) for paths in path_list for path in paths]
        weight_avg = [path_average_weight(graph, path) for paths in path_list for path in paths]
        graph = select_best_path(graph, [path for paths in path_list for path in paths], path_length, weight_avg, delete_entry_node=True, delete_sink_node=False)
    return graph


def solve_out_tips(graph, ending_nodes):
    '''
    Selects the best sink for the graph, based on all the possible sinks.
    '''
    path_list = [list(nx.all_simple_paths(graph, node, end)) for node in graph.nodes for end in ending_nodes if len(list(graph.successors(node))) > 1]
    if path_list:
        path_length = [len(path) for paths in path_list for path in paths]
        weight_avg = [path_average_weight(graph, path) for paths in path_list for path in paths]
        graph = select_best_path(graph, [path for paths in path_list for path in paths], path_length, weight_avg, delete_entry_node=False, delete_sink_node=True)
    return graph


def get_starting_nodes(graph):
    """Get nodes without predecessors in the graph.

    :param graph: (nx.DiGraph) A directed graph object.
    :return: A list of nodes without predecessors.
    """
    starting_nodes = [node for node in graph.nodes() if not any(graph.predecessors(node))]
    return starting_nodes

def get_sink_nodes(graph):
    """Get nodes without successors in the graph.

    :param graph: (nx.DiGraph) A directed graph object.
    :return: A list of nodes without successors.
    """
    sink_nodes = [node for node in graph.nodes() if not any(graph.successors(node))]
    return sink_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract contigs from the graph.

    :param graph: (nx.DiGraph) A directed graph object.
    :param starting_nodes: A list of nodes without predecessors.
    :param ending_nodes: A list of nodes without successors.
    :return: A list of tuples (contig, contig_length).
    """
    contig_list=[]
    
    for start_node in starting_nodes:
        for end_node in ending_nodes:
                for path in nx.all_simple_paths(graph, source=start_node, target=end_node):
                     contig=path[0]

                for node in path[1:]:
                    contig=contig+node[-1]

                contig_size=len(contig)
                contig_list.append((contig, contig_size))

    return contig_list       


def save_contigs(contigs_list, output_file):
    """Write contigs to a FASTA file.

    :param contigs_list: A list of tuples (contig, contig_length).
    :param output_file: (str) Path to the output file.
    """
    with open(output_file, "w") as f:
        for i, (contig, length) in enumerate(contigs_list):
            header = f">contig_{i} len={length}\n"
            wrapped_contig = textwrap.fill(contig, width=80)
            f.write(header + wrapped_contig + "\n")


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file.resolve())


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """

     # Get arguments
    args = get_arguments()
    fastq=args.fastq_file
    kmer_length=args.kmer_size
    output=os.curdir + os.sep+args.output_file

    # Read the file and build a graph
    print('building a graph...')
    kmer_dict=build_kmer_dict(fastq,kmer_length)
    graph=build_graph(kmer_dict)
    print (graph)

    # Resolving the bubbles
    print('resolving the bubbles...')
    graph=simplify_bubbles(graph)

    #Resolving the starting and sink nodes
    print('resolving the starting and sink nodes...')
    print('starting nodes...')
    starting_nodes=get_starting_nodes(graph)
    print('sink nodes...')
    sink_nodes=get_sink_nodes(graph)
    print('solving entry...')
    graph=solve_entry_tips(graph,starting_nodes)
    print('solving out...')
    graph=solve_out_tips(graph,sink_nodes)

    #Writing the contig(s)
    print('saving the graph to png...')

    labels = nx.get_edge_attributes(graph,'weight')
    positions=nx.spring_layout(graph, scale=8)
    nx.draw(graph, pos=positions, node_size=10)
    nx.draw_networkx_edge_labels(graph, pos=positions, edge_labels=labels, font_size=6)
    plt.savefig('Graph.png')

    print('saving the contigs to fasta...')
    starting_nodes=get_starting_nodes(graph)
    sink_nodes=get_sink_nodes(graph)
    contig_list=get_contigs(graph, starting_nodes, sink_nodes)
    output_file = "contigs.fasta"  
    save_contigs(contig_list,output)

if __name__ == '__main__':
    main()








