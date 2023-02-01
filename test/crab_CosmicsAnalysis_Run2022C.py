import CRABClient
from CRABClient.UserUtilities import config

config = config()

# General
config.General.workArea = '/eos/user/r/rlopezru/Run2022C-PromptReco-v1/MINIAOD/crab_projects'
config.General.requestName = 'CosmicsAnalysis_Run2022C'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runNtuplizer_cfg.py'
config.JobType.maxMemoryMB = 2500
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['Cosmics_Run2022C.root']

# Data
config.Data.inputDataset = '/NoBPTX/Run2022C-PromptReco-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100 # splitting to get ~500 output files
config.Data.outLFNDirBase = '/store/user/rlopezru/Cosmics' # modify accordingly
config.Data.publication = False
config.Data.outputDatasetTag = 'CosmicsAnalysis_Run2022C'

# Site
config.Site.storageSite = 'T3_CH_CERNBOX'
