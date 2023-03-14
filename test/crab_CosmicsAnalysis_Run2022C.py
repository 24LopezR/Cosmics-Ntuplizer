import CRABClient
from CRABClient.UserUtilities import config, getLumiListInValidFiles
#from FWCore.DataStructs.LumiList import LumiList
from FWCore.PythonUtilities.LumiList import LumiList

config = config()

# General
config.General.workArea = '/eos/user/r/rlopezru/crab_projects'
config.General.requestName = 'CosmicsAnalysis_Run2022C_AOD-Ntuples'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Cosmics_runNtuplizer_cfg.py'
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['Cosmics_Run2022C_AOD-Ntuples.root']

# Data
config.Data.inputDataset = '/NoBPTX/Run2022C-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100 # splitting to get ~500 output files
config.Data.outLFNDirBase = '/store/user/rlopezru/Samples/' # modify accordingly
config.Data.publication = False
config.Data.outputDatasetTag = 'CosmicsAnalysis_Run2022C_AOD-Ntuples'

# LumiList
'''
json_dir = '/eos/user/c/cmsdqm/www/CAF/certification/Cosmics22/'
json_list = ['Cosmics22_CRAFT_349078_349529_call4_GoldenJSON.txt',
             'Cosmics22_CRAFT_349610_349839_call5_GoldenJSON.txt',
             'Cosmics22_CRAFT_349840_350164_call6_GoldenJSON.txt',
             'Cosmics22_CRAFT_350166_350457_call7_GoldenJSON.txt']
json_list = [json_dir+f for f in json_list]
print(json_list)
lumi_list = [LumiList(filename=json_dir+'Cosmics22_CRAFT_349078_349529_call4_GoldenJSON.txt'),
             LumiList(filename=json_dir+'Cosmics22_CRAFT_349610_349839_call5_GoldenJSON.txt')]
full_lumi_list = lumi_list[0] + lumi_list[1]
''' 
#config.Data.lumiMask = '/afs/cern.ch/user/r/rlopezru/private/ntuplizer_test/CMSSW_12_4_0/src/Analysis/Cosmics-Ntuplizer/full_lumi_mask.json'

# Site
config.Site.storageSite = 'T3_CH_CERNBOX'
