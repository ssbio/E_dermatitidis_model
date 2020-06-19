#! /usr/bin/python
#! python3.7
#try to specify that will use python version 3.7
__author__      = "Wheaton Schroeder"

#This code is created to read the gams-based input files for a genome scale model and instead create a GSM using the 
#cobrapy tool. This is then followed by performing flux balance analysis (FBA) on the given model

#This was written to satisfy the request of reviewers for the protocol paper:
#Optfill, Shadow Price, and the Genome-Scale Exophiala dermatitidis Reconstruction iEde2091: Procedure and Application
#as a Model of Pigment Synthesis
#By: Wheaton L Schroeder and Rajib Saha
#Submitted to STAR Protocols

#latest version: 06/17/2020

#run the model convert files
import os
import cobra

#from cobra import Model, Reaction, Metabolite
from cobra import Model, Reaction, Metabolite

#read the model with cobrapy
model = cobra.io.read_sbml_model("iEde2091.xml")

#output file, report things like reaction rates
output=open('FBA_python_output.tsv','w')

#specify the bounds
M = 10000       #positive really big number (arbitrary)

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = -0.25      #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = 0          #glucose
model.reactions.get_by_id("R012ex").lower_bound = 0          #urate

model.reactions.get_by_id("R001ex").upper_bound = M          #water
model.reactions.get_by_id("R002ex").upper_bound = 0          #phosphate
model.reactions.get_by_id("R003ex").upper_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").upper_bound = 0          #ammonia
model.reactions.get_by_id("R005ex").upper_bound = 0          #sulfate
model.reactions.get_by_id("R006ex").upper_bound = M          #protons
model.reactions.get_by_id("R007ex").upper_bound = M          #carbon dioxide
model.reactions.get_by_id("R008ex").upper_bound = 0          #oxygen
model.reactions.get_by_id("R009ex").upper_bound = M          #ethanol
model.reactions.get_by_id("R010ex").upper_bound = 0          #acetate
model.reactions.get_by_id("R011ex").upper_bound = 0          #glucose
model.reactions.get_by_id("R012ex").upper_bound = M          #urate

model.repair()

solution = model.optimize()

print("\n\nobjective value:    "+str(solution.objective_value)+" h^-1\n\n")
print("water exchange:     "+str(solution.fluxes["R001ex"])+" h^-1")
print("phosphate exchange: "+str(solution.fluxes["R002ex"])+" h^-1")
print("sucrose exchange:   "+str(solution.fluxes["R003ex"])+" h^-1")
print("ammonia exchange:   "+str(solution.fluxes["R004ex"])+" h^-1")
print("sulfate exchange:   "+str(solution.fluxes["R005ex"])+" h^-1")
print("proton exchange:    "+str(solution.fluxes["R006ex"])+" h^-1")
print("CO2 exchange:       "+str(solution.fluxes["R007ex"])+" h^-1")
print("oxygen exchange:    "+str(solution.fluxes["R008ex"])+" h^-1")
print("enthanol exchange:  "+str(solution.fluxes["R009ex"])+" h^-1")
print("acetate exchange:   "+str(solution.fluxes["R010ex"])+" h^-1")
print("glucose exchange:   "+str(solution.fluxes["R011ex"])+" h^-1")
print("urate exchange:     "+str(solution.fluxes["R012ex"])+" h^-1\n")

output.write("reaction"+"\t"+"FBA level\n")

for rxn in model.reactions:	        #for each line in the input file    

    output.write(rxn.id+"\t"+str(solution.fluxes[rxn.id])+"\n")

