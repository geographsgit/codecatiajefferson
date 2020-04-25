#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 19:40:24 2020

@author: catia(souz.kti@gmail.com) and Jefferson
"""

from igraph import Graph
#import tkinter as tk
#from tkinter import filedialog
import numpy as np
import pandas as pd
from statistics import mean 
import math


class GraphMaking:
    
    COLUMNS = [
        'limiar',       
        'avg_betweenness_eccentricity',
        'avg_path_length',
        'avg_degree',
        'diameter',
        'transitivity',
        'number_of_edges'
    ]
    
    COLUMNS_NODE = ['degree', 'betweenness','transitivity', 'shortest_paths', 
                    'limiar', 'diameter_graph', 'number_of_edges_graph',
                    'avg_degree','avg_betweenness_eccentricity',
                    'transitivity','avg_path_length']

    def __init__(self):
        self.graph = None
        self.flow = []
        self.contador = 0
        self.order = 0
        self.infinities = []
        self.vuln = []
        self.weighted_vuln = []

        self.pop = 0
        self.ek = []
        self.infinities_weight = []

        for line in self.flow:
            self.pop += sum(line)

    def create_graph(self, filename):
        '''Creates a graph from adjacency matrix'''
        
        df=pd.read_csv(filename, sep=',',header=None, skiprows=1)
        df = df.drop(df.columns[0], axis=1)
        
        A=np.matrix(df.values)
        self.graph = Graph.Adjacency(A.tolist())
        self.graph = self.graph.as_undirected()

        self.order = self.graph.vcount()
        self.n = self.order

    def global_efficiency(self, g):
        invcam = 0
        self.mencam = g.shortest_paths_dijkstra()  # Função que retorna o menor caminho entre cada par de vértices do grafo
        # 'If' com o objetivo de manter 'n' constante (pois após remover um vértice n se torna n-1)
        for vertice in self.mencam:  # itera para cada vertice de origem na lista de menores caminhos
            for caminho in vertice:  # itera para cada destino do vertice de origem em questão
                if caminho != 0:
                    invcam += 1 / caminho  # Acumula os valores de eficiencia associada a cada par de vértices
        eg = invcam / (self.n * (self.n - 1))
        if self.contador == 0:
            self.origi_efi = eg
        else:
            self.ek.append(eg)
        self.contador += 1
        return eg

    def weighted_global_efficiency(self, g):
        invcam = 0
        pop = 0
        for lista in self.flow:
            for elemento in lista:
                pop += elemento
        self.mencam = g.shortest_paths_dijkstra()  # Função que retorna o menor caminho entre cada par de vértices do grafo
        # 'If' com o objetivo de manter 'n' constante (pois após remover um vértice n se torna n-1)
        for i, vertice in enumerate(self.mencam):  # itera para cada vertice de origem na lista de menores caminhos
            for j, caminho in enumerate(vertice):  # itera para cada destino do vertice de origem em questão
                if caminho != 0:
                    invcam += (1 / caminho) * self.flow[i][j]  # Acumula os valores de eficiencia associada a cada par de vértices

        eg = invcam / (self.n * (self.n - 1) * pop)
        return eg

    def vulnerability(self):
        # Eficiencia com o vertice
        eg = self.global_efficiency(self.graph)
        # Eficiencia sem o vertice
        for i in range(0, self.order):
            g = self.graph.copy()  # A cada iteração a varoável 'graf' receberá a cópia do grafo original
            del_list = []  # lista que conterá os ids das arestas do vértice 'i' que serão removidas.
            for target_vertex_id in range(0, self.order):
                try:
                    del_list.append(g.get_eid(i,target_vertex_id))  # função que captura o id da aresta pertencente ao par de vértices (i, target_vertex_id), e acrescenta na lista 'del_list'
                except:
                    pass  # caso o id não exista
            g.delete_edges(del_list)
            efi = self.global_efficiency(g)
            v = (eg - efi) / eg
            self.vuln.append(v)

    def weighted_vulnerability(self):
        # Eficiencia com o vertice
        eg = self.weighted_global_efficiency(self.graph)
        # Eficiencia sem o vertice
        for i in range(0, self.order):
            g = self.graph.copy()  # A cada iteração a varoável 'graf' receberá a cópia do grafo original
            del_list = []  # lista que conterá os ids das arestas do vértice 'i' que serão removidas.
            for target_vertex_id in range(0, self.order):
                try:
                    del_list.append(g.get_eid(i, target_vertex_id))  # função que captura o id da aresta pertencente ao par de vértices (i, target_vertex_id), e acrescenta na lista 'del_list'
                except:
                    pass  # caso o id não exista
            g.delete_edges(del_list)
            efi = self.weighted_global_efficiency(g)
            v_w = (eg - efi) / eg
            self.weighted_vuln.append(v_w)

    def isolation(self):
        for i in range(0,self.order):  # for de 0 até a quantidade de vértices, onde cada iteração representa uma rua que ficará inacessível a todas as outras.
            self.acumula_infinito = 0  # Variável que guardará a quantidade de infinitos após as arestas do vértice 'i' ser removida
            graph_copy = self.graph.copy()  # A cada iteração a varoável 'graf' receberá a cópia do grafo original
            del_list = []  # lista que conterá os ids das arestas do vértice 'i' que serão removidas.
            for target_vertex_id in range(0, self.order):
                try:
                    del_list.append(graph_copy.get_eid(i,target_vertex_id))  # função que captura o id da aresta pertencente ao par de vértices (i, target_vertex_id), e acrescenta na lista 'del_list'
                except:
                    pass  # caso o id não exista
            graph_copy.delete_edges(del_list)  # deleta as arestas contidas na lista 'del_list'
            self.mencam = graph_copy.shortest_paths_dijkstra()  # cria lista com todos os menores caminhos para de/para todos os vértices
            for vertice in self.mencam:  # itera para cada vertice de origem na lista de menores caminhos
                self.acumula_infinito += vertice.count(float('inf'))  # Guarda as ruas inacessiveis para um determinado vértice quando ele é removido
            self.infinities.append(self.acumula_infinito)  # lista que guarda em cada posição a quantidade de infinitos relacionados ao isolamento da rua 'i'.


    def weighted_isolation(self):
        for i in range(0,self.order):  # for de 0 até a quantidade de vértices, onde cada iteração representa uma rua que ficará inacessível a todas as outras.
            self.acumula_infinito = 0  # Variável que guardará a quantidade de infinitos após as arestas do vértice 'i' ser removida
            graph_copy = self.graph.copy()  # A cada iteração a varoável 'graf' receberá a cópia do grafo original
            del_list = []  # lista que conterá os ids das arestas do vértice 'i' que serão removidas.
            for target_vertex_id in range(0, self.order):
                try:
                    del_list.append(graph_copy.get_eid(i,target_vertex_id))  # função que captura o id da aresta pertencente ao par de vértices (i, target_vertex_id), e acrescenta na lista 'del_list'
                except:
                    pass  # caso o id não exista
            graph_copy.delete_edges(del_list)  # deleta as arestas contidas na lista 'del_list'
            self.mencam = graph_copy.shortest_paths_dijkstra()  # cria lista com todos os menores caminhos para de/para todos os vértices
            for i, vertice in enumerate(self.mencam):  # itera para cada vertice de origem na lista de menores caminhos
                for j, path in enumerate(vertice):
                    if path==float('inf'):
                        self.acumula_infinito += (1*self.flow[i][j])/self.pop
            self.infinities_weight.append(self.acumula_infinito)  # lista que guarda em cada posição a quantidade de infinitos relacionados ao isolamento da rua 'i'.

    def create_results(self, grp):
        output_global = open('output_global.txt', 'w')
        output_local = open('output_local.txt', 'a')
        output_global.write('N, pop, Original_Efi, \n' + str(self.order) + ',' + str(self.pop) + ',' + str(self.origi_efi))
        output_local.write('id,Ek,Vul_efi,Vul_efi_w,vul_iso,vul_iso_w \n')
        for i in range(0, grp.order):
            output_local.write(str(self.vs[i]) + ',' + str(self.ek[i]) + ',' + str(self.vuln[i]) + ',' + str(self.weighted_vuln[i]) + ',' + str(self.infinities[i])+ ','+str(self.infinities_weight[i]) +'\n')
   
    def betweenness_eccentricity(self):
        '''Returns avg_betweenness_eccentricity'''
        return sum(self.graph.betweenness(vertices=None, directed=False, cutoff=None, weights=None, nobigint=True))/(self.graph.vcount())
    
    def metrics_of_graph(self, limiar):
        '''Returns graph metrics'''
        betweenness_eccentricity = self.betweenness_eccentricity()
        average_path_length = self.graph.average_path_length(directed=False, unconn=True)
        average_degree = mean(self.graph.degree())
        diameter = self.graph.diameter(directed=False, unconn=True, weights=None)
        transitivity = self.graph.transitivity_undirected(mode="nan")
        number_of_edges = self.graph.ecount()
        
       
        metrics = pd.DataFrame(columns=self.COLUMNS)
        list_of_metrics = [limiar, betweenness_eccentricity,
            average_path_length,
            average_degree, 
            diameter,
            transitivity,
            number_of_edges],
                
        row_metrics = pd.DataFrame(list_of_metrics,columns=self.COLUMNS)
        
        metrics = pd.concat([metrics, row_metrics])
        metrics.to_csv('metrics_'+str(limiar)+'.csv')
       
    def metrics_of_nodes(self, limiar):
        '''Returns metrics for all nodes'''
        
        number_of_nodes = len(self.graph.vs)
        
        shortest_paths=[0]*number_of_nodes
        
        pos=0
        for lista in self.graph.shortest_paths_dijkstra():
            count=0
            for elem in lista:
                if(np.isfinite(elem)):
                    count+=elem
            shortest_paths[pos]= (count/(number_of_nodes-1))#Conta o proprio elemento?
            pos+=1
        
        betweenness_eccentricity = self.betweenness_eccentricity()
        average_path_length = self.graph.average_path_length(directed=False, unconn=True)
        average_degree = mean(self.graph.degree())
        diameter = self.graph.diameter(directed=False, unconn=True, weights=None)
        transitivity = self.graph.transitivity_undirected(mode="nan")
        number_of_edges = self.graph.ecount()
        
        transitivities = [-1 if math.isnan(x) else x for x in self.graph.transitivity_local_undirected(mode='nan')]
        
        df = pd.DataFrame(self.graph.degree())
        df2 = pd.DataFrame(self.graph.betweenness(vertices=None, directed=False, cutoff=None, weights=None, nobigint=True))
        df3 = pd.DataFrame(transitivities)
        df4 = pd.DataFrame(shortest_paths)
        df5 = pd.DataFrame([limiar] * number_of_nodes)
        df6 = pd.DataFrame([diameter] * number_of_nodes)
        df7 = pd.DataFrame([number_of_edges] * number_of_nodes)
        df8 = pd.DataFrame([average_degree] * number_of_nodes)
        df9 = pd.DataFrame([betweenness_eccentricity]* number_of_nodes)
        df10 = pd.DataFrame([transitivity] * number_of_nodes)
        df11 = pd.DataFrame([average_path_length] * number_of_nodes)
        
        node_metrics = pd.concat([df, df2], axis=1, sort=False)
        node_metrics = pd.concat([node_metrics, df3], axis=1, sort=False)
        node_metrics = pd.concat([node_metrics, df4], axis=1, sort=False)
        node_metrics = pd.concat([node_metrics, df5], axis=1, sort=False)
        node_metrics = pd.concat([node_metrics, df6], axis=1, sort=False)
        node_metrics = pd.concat([node_metrics, df7], axis=1, sort=False)
        node_metrics = pd.concat([node_metrics, df8], axis=1, sort=False)
        node_metrics = pd.concat([node_metrics, df9], axis=1, sort=False)
        node_metrics = pd.concat([node_metrics, df10], axis=1, sort=False)
        node_metrics = pd.concat([node_metrics, df11], axis=1, sort=False)
        
        
        node_metrics.columns = self.COLUMNS_NODE

        node_metrics.to_csv('node_metrics'+str(limiar)+'.csv')

    def main(self):
        print (self.create_graph.__doc__)
        limiar = 0
#        filename = 'adj_matrix_roads'+str(limiar)+'.csv'
        filename = 'simetric_matrix.csv'

        self.create_graph(filename)
        self.metrics_of_graph(limiar)
        self.metrics_of_nodes(limiar)

        self.create_results(self)
        self.vulnerability()
        self.isolation()
        self.weighted_isolation()
        self.weighted_vulnerability()
        
if __name__ == '__main__':
    grp = GraphMaking()
    grp.main()
        