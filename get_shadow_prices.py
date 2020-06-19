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
output=open('shadow_price_python_output.tsv','w')

#need to write the header for shadow costs
output.write("Trial#\tcode\tgrowth\tX00001[c]\tC04033[c]\tC00779[c]\tC01173[c]\tC01624[c]\tC00083[c]\tC17937DHN[e]\t")
output.write("C05578[e]\tC05579[e]\tC01693[e]\tC05604[e]\tC00355[c]\tC00822[c]\tC00082[c]\tC17937Eu[e]\tC00544[c]\t")
output.write("X00002[c]\tC01179[c]\tC17937Pyo[e]\tC00353[c]\tC03427[c]\tC05421[c]\tC05414[c]\tC05430[c]\tC05431[c]\t")
output.write("C05432[c]\tC05435[c]\tC15867[c]\tC08613[c]\tC02094[c]\tC19892[c]\tC08607[c]\tC00097[c]\tModel Status\n")

#specify the bounds
M = 10000       #positive really big number (arbitrary)

#track trial number
trial = 1;

########################################## TRIAL 1 ##########################################
#            Carbon-limited sucrose as carbon source, LB(R003ex) = -0.25                    #
#############################################################################################

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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("SCL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 2 ##########################################
#            Carbon-limited sucrose as carbon source, LB(R003ex) = -0.50                    #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = -0.50      #sucrose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("SCM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 3 ##########################################
#            Carbon-limited sucrose as carbon source, LB(R003ex) = -1                       #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = -1         #sucrose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("SCH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 4 ##########################################
#            Carbon-limited ethanol as carbon source, LB(R009ex) = -1.5                     #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = -1.5       #ethanol
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ECL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 5 ##########################################
#            Carbon-limited ethanol as carbon source, LB(R009ex) = -3                       #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = -3         #ethanol
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ECM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 6 ##########################################
#            Carbon-limited ethanol as carbon source, LB(R009ex) = -6                       #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = -6         #ethanol
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ECH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 7 ##########################################
#            Carbon-limited acetate as carbon source, LB(R010ex) = -1.5                     #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = -1.5       #acetate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ACL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 8 ##########################################
#            Carbon-limited acetate as carbon source, LB(R010ex) = -3                      #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = -3         #acetate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ACM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 9 ##########################################
#            Carbon-limited acetate as carbon source, LB(R010ex) = -6                      #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = -6         #acetate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ACH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 10 #########################################
#            Carbon-limited glucose as carbon source, LB(R011ex) = -0.5                     #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = -0.5       #glucose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("GCL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 11 #########################################
#            Carbon-limited glucose as carbon source, LB(R011ex) = -1                       #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = -1         #glucose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("GCM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 12 #########################################
#            Carbon-limited glucose as carbon source, LB(R011ex) = -2                       #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = -2         #glucose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("GCH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial

########################################## TRIAL 13 #########################################
#            Nitrogen-limited sucrose as carbon source, LB(R004ex) = -0.05                  #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = -M/12       #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.05      #ammonia
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("SNL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 14 #########################################
#            Nitrogen-limited sucrose as carbon source, LB(R004ex) = -0.1                   #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = -M/12       #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.1       #ammonia
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("SNM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 15 #########################################
#            Nitrogen-limited sucrose as carbon source, LB(R004ex) = -0.2                   #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = -M/12      #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.2       #ammonia
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("SNH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 16 #########################################
#            Nitrogen-limited ethanol as carbon source, LB(R004ex) = -0.05                  #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.05      #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = -M/2        #ethanol
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ENL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 17 #########################################
#            Nitrogen-limited ethanol as carbon source, LB(R004ex) = -0.1                   #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.1       #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = -M/2        #ethanol
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ENM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 18 #########################################
#            Nitrogen-limited ethanol as carbon source, LB(R004ex) = -0.2                   #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.2       #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = -M/2        #ethanol
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ENH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 19 #########################################
#            Nitrogen-limited acetate as carbon source, LB(R004ex) = -0.05                  #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.05      #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = -M/2        #acetate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ANL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 20 #########################################
#            Nitrogen-limited acetate as carbon source, LB(R004ex) = -0.1                   #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.1       #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = -M/2        #acetate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ANM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 21 #########################################
#            Nitrogen-limited acetate as carbon source, LB(R004ex) = -0.2                   #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.2       #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = -M/2        #acetate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ANH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 22 #########################################
#            Nitrogen-limited glucose as carbon source, LB(R004ex) = -0.05                  #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.05      #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = -M/6        #glucose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("GNL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 23 #########################################
#            Nitrogen-limited glucose as carbon source, LB(R004ex) = -0.1                   #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.1       #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = -M/6        #glucose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("GNM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 23 #########################################
#              Sulfur-limited glucose as carbon source, LB(R004ex) = -0.2                   #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -0.1       #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -M         #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = 0          #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = -M/6        #glucose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("GNH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 25 #########################################
#              Sulfur-limited glucose as carbon source, LB(R005ex) = -0.0005                #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = -M/12      #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.0005    #sulfate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("SSL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 26 #########################################
#              Sulfur-limited glucose as carbon source, LB(R005ex) = -0.001                 #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = -M/12      #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.001    #sulfate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("SSM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 27 #########################################
#              Sulfur-limited glucose as carbon source, LB(R005ex) = -0.002                 #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = -M/12      #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.002     #sulfate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("SSH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 28 #########################################
#              Sulfur-limited ethanol as carbon source, LB(R005ex) = -0.0005                #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.0005    #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = -M/2       #ethanol
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ESL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 29 #########################################
#              Sulfur-limited ethanol as carbon source, LB(R005ex) = -0.001                 #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.001     #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = -M/2       #ethanol
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ESM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 30 #########################################
#              Sulfur-limited ethanol as carbon source, LB(R005ex) = -0.002                 #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.002     #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound = -M/2       #ethanol
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ESH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 31 #########################################
#              Sulfur-limited acetate as carbon source, LB(R005ex) = -0.0005                #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.0005    #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound =  0         #ethanol
model.reactions.get_by_id("R010ex").lower_bound = -M/2       #acetate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ASL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 32 #########################################
#              Sulfur-limited acetate as carbon source, LB(R005ex) = -0.001                 #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.001     #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound =  0         #ethanol
model.reactions.get_by_id("R010ex").lower_bound = -M/2       #acetate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ASM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 33 #########################################
#              Sulfur-limited acetate as carbon source, LB(R005ex) = -0.002                 #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.002     #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound =  0         #ethanol
model.reactions.get_by_id("R010ex").lower_bound = -M/2       #acetate
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("ASH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 34 #########################################
#              Sulfur-limited glucose as carbon source, LB(R005ex) = -0.0005                #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.0005    #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound =  0         #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = -M/6       #glucose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("GSL\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 35 #########################################
#              Sulfur-limited glucose as carbon source, LB(R005ex) = -0.001                 #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.001     #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound =  0         #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = -M/6       #glucose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("GSM\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

########################################## TRIAL 36 #########################################
#              Sulfur-limited glucose as carbon source, LB(R005ex) = -0.002                 #
#############################################################################################

#set the bounds for exchange reactions
#note that ".lowerbound" does not appear to work, only ".lower_bound"
#note that ".upperbound" does not appear to work, only ".upper_bound"
model.reactions.get_by_id("R001ex").lower_bound = -M         #water
model.reactions.get_by_id("R002ex").lower_bound = -M         #phosphate
model.reactions.get_by_id("R003ex").lower_bound = 0          #sucrose
model.reactions.get_by_id("R004ex").lower_bound = -M         #ammonia
model.reactions.get_by_id("R005ex").lower_bound = -0.002     #sulfate
model.reactions.get_by_id("R006ex").lower_bound = -M         #protons
model.reactions.get_by_id("R007ex").lower_bound = 0          #carbon dioxide
model.reactions.get_by_id("R008ex").lower_bound = -M         #oxygen
model.reactions.get_by_id("R009ex").lower_bound =  0         #ethanol
model.reactions.get_by_id("R010ex").lower_bound = 0          #acetate
model.reactions.get_by_id("R011ex").lower_bound = -M/6       #glucose
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

#solve FBA
solution = model.optimize()

#need to write the shadow prices
output.write(str(trial)+"\t")
output.write(str("GSH\t"))
output.write("{:e}".format(solution.objective_value)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00001__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C04033__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00779__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01173__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01624__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00083__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937DHN__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05578__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05579__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01693__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05604__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00355__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00822__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00082__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Eu__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00544__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.X00002__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C01179__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C17937Pyo__e)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00353__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C03427__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05421__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05414__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05430__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05431__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05432__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C05435__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C15867__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08613__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C02094__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C19892__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C08607__c)+"\t")
output.write("{:e}".format(solution.shadow_prices.C00097__c)+"\t")
output.write("1\n")

trial=trial+1
print("Trial "+str(trial)+" running...\n")

print("SUCCESS")