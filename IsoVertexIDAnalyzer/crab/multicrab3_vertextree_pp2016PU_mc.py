if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    from CRABClient.UserUtilities import config, getUsernameFromSiteDB
    config = config()

    config.General.workArea = 'Ana'
    config.General.transferOutputs = True
    config.General.transferLogs = False
    config.JobType.pluginName = 'Analysis'
    config.JobType.allowUndistributedCMSSW = True
#    config.JobType.maxMemoryMB = 3000
#    config.Data.unitsPerJob = 10
    config.Data.totalUnits = 200
#    config.Data.inputDBS = 'phys03'
    config.Data.splitting = 'FileBased'
#    config.Data.lumiMask = 'json_pPb.txt'
#    config.Data.lumiMask = 'Cert_285479-286496_HI8TeV_PromptReco_pPbPbp_Collisions16_JSON_NoL1T.txt'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.publication = False
#    config.Site.storageSite = 'T2_US_MIT'
    config.Site.storageSite = 'T3_US_Rice'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = 'vertex_MC16NoPU_v2'
    config.JobType.psetName = '../test/vertextree_mc_cfg.py'
    config.Data.inputDataset = '/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIISummer16DR80-NoPU_80X_mcRun2_asymptotic_v14-v1/AODSIM'
    config.Data.unitsPerJob = 1
    submit(config)
