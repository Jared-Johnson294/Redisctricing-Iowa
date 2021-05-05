#!/usr/bin/env python
# coding: utf-8

# In[1]:


from gerrychain import Graph


# In[2]:


#read the iowa json county graph
filepath = 'C:\\district-data\\'
filename = 'COUNTY_19.json'
G = Graph.from_json( filepath + filename )


# In[3]:


#Make sure the file is read correctly by printing the node #, and its population
for node in G.nodes:
    name = G.nodes[node]["NAME10"]
    population = G.nodes[node]['TOTPOP']
    x_coordinate = G.nodes[node]['C_X']
    y_coordinate = G.nodes[node]['C_Y']
    print("Node",node,"is",name,"County, which has population on",population,"and is centered at (",x_coordinate,",",y_coordinate,")")


# In[4]:


#create distance dictionary
from geopy.distance import geodesic

Marshal = ( G.nodes[0]['C_Y'], G.nodes[0]['C_X'] )
Boone = ( G.nodes[52]['C_Y'], G.nodes[52]['C_X'] )
Polk = ( G.nodes[73]['C_Y'], G.nodes[73]['C_X'] )

print("Marshal -> Boone",geodesic(Marshal, Boone).miles)
print("Boone -> Polk", geodesic(Boone, Polk).miles)
print("Polk -> Marshal",geodesic(Polk, Marshal).miles)


# In[5]:


#Check the dictionary 
dist = dict()
for i in G.nodes:
    for j in G.nodes:
        loc_i = ( G.nodes[i]['C_Y'], G.nodes[i]['C_X'] )
        loc_j = ( G.nodes[j]['C_Y'], G.nodes[j]['C_X'] )
        dist[i,j] = geodesic(loc_i,loc_j).miles


# In[6]:


print("Marshal -> Boone:",dist[0,52])


# In[7]:


#Find the upper and lower bound 
deviation = 0.01

import math 
K = 4
total_population = sum(G.nodes[node]['TOTPOP'] for node in G.nodes)

L = math.ceil((1-deviation/2)*total_population/K)
U = math.floor((1+deviation/2)*total_population/K)
print("Using L =",L,"and U =",U,"and K =",K)


# In[8]:


#Start of the model
import gurobipy as gp
from gurobipy import GRB

m = gp.Model()
#X[i,j] variable which equals one when county i 
#       Is assigned to (the district centered at) county j
X = m.addVars(G.nodes, G.nodes, vtype=GRB.BINARY)


# In[9]:


#Objective equation
m.setObjective( gp.quicksum( dist[i,j]*dist[i,j]*G.nodes[i]['TOTPOP']*X[i,j] for i in G.nodes for j in G.nodes), GRB.MINIMIZE )


# In[10]:


#Constraint equations
#each county i is assigned to one district
m.addConstrs( gp.quicksum(X[i,j] for j in G.nodes) == 1 for i in G.nodes) 
#there should be K 
m.addConstr( gp.quicksum( X[j,j] for j in G.nodes ) == K )
#the population in each district has to be bwtween the upper and lower limits
m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * X[i,j] for i in G.nodes) >= L * X[j,j] for j in G.nodes )
m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * X[i,j] for i in G.nodes) <= U * X[j,j] for j in G.nodes )
#if i is assigned to j then j is center
m.addConstrs( X[i,j] <= X[j,j] for i in G.nodes for j in G.nodes )

m.update()


# In[11]:


#Add the contiguity by adding mpre constraints
import networkx as nx 
DG = nx.DiGraph(G)
#Define new variables 
f = m.addVars( DG.nodes, DG.edges, vtype=GRB.CONTINUOUS)
M = DG.number_of_nodes()-1
#node j cannot recieve flow of its own type
m.addConstrs( gp.quicksum( f[j,u,j] for u in DG.neighbors(j) ) == 0 for j in DG.nodes )
#node i can recieve flow of type j only if i is assigned to j
m.addConstrs( gp.quicksum( f[j,u,i] for u in DG.neighbors(i)) <= M * X[i,j] for i in DG.nodes for j in DG.nodes if i != j )
#i is assigned to j, then i should consume one unit of j flow 
#     Otherwise,, i should consume no units 
m.addConstrs( gp.quicksum( f[j,u,i] - f[j,i,u] for u in DG.neighbors(i)) == X[i,j] for i in DG.nodes for j in DG.nodes if i != j )

m.update()


# In[12]:


#optimize
m.Params.MIPGap = 0.0
m.optimize()


# In[13]:


#Print what counties are in each district
print("The moment of inertia objective is",m.objval)

centers = [j for j in G.nodes if X[j,j].x > 0.5 ]
districts = [ [i for i in G.nodes if X[i,j].x > 0.5] for j in centers]
district_counties = [ [ G.nodes[i]["NAME10"] for i in districts[j] ] for j in range(K)]
district_populations = [ sum(G.nodes[i]["TOTPOP"] for i in districts[j]) for j in range(K) ]
for j in range(K):
    print("District",j,"has population",district_populations[j],"and contains counties",district_counties[j])


# In[14]:


#Start creation of the map
import geopandas as gpd


# In[15]:


#Find files in computer to graph iowa
filepath_ = 'C:\\district-data\\'
filename_ = 'IA_counties.shp'

df = gpd.read_file( filepath_ + filename_ )


# In[16]:


#Finish maps
assignment = [ -1 for u in G.nodes ]

for j in range(len(districts)):
    
    for i in districts[j]:
    
        geoID = G.nodes[i]["GEOID10"]
    
        for u in G.nodes:
            if geoID == df["GEOID10"][u]:
                assignment[u] = j
        
df['assignment'] = assignment
my_fig = df.plot(column='assignment').get_figure()

