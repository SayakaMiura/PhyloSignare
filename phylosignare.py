from Mutation_profile.trinucleotide import trinucleotide
from parsers.ParamsParser import ParamsParser
from decomposition.Exposure_computer import Exposure_computer
from decomposition.significanse_test import significanse_test
from parsers.ParamsParserSig import ParamsParserSig
from output.file_writer import file_writer
import glob
import os
import sys
import shutil

Input=sys.argv[1] #*.inpu
Control=sys.argv[2]
ChiTa='chisquareDiffTestTa.txt'
print('reading control file (Control.txt)')
ParamSig=ParamsParserSig(Control)
COS,Rpath,Method,Can2sig,Expected_SigLsLs,EstimateExposureRtemplate=ParamSig.get_parameter_setting()
summary_out=ParamSig.get_summary_out()
print(summary_out)
print(COS,Rpath,Method,Can2sig,Expected_SigLsLs,EstimateExposureRtemplate)
MinMutC=1 #20: make combi for TRACEx
OutPutCut=20
Fol=Input.replace(Input.split('\\')[-1],'')
SigDif=0.02
BigSig=0.4
CutSig=0.01
scoreCut=0.5
FigOut='n'
if len(sys.argv)>3: FigOut='y'
Expected_SigLs=','.join(map(str,Expected_SigLsLs))
OutFolName=list(Can2sig.keys())[0]+'-PhyloSignare' 
COSin=open(COS,'r').readlines()[0].split(',')
TotSig=len(COSin)-1
print(TotSig)
cwd = os.getcwd()  
os.chdir(Fol)
Res=Input.split('\\')[-1][:-6]+'-'+OutFolName
if os.path.isdir(Res)!=True: os.mkdir(Res) 
os.chdir(cwd)	
ResDir=Fol+Res+'\\'
print('1. combine or filter short branches')
Param=ParamsParser(Input)
SaveFile=file_writer()
Format_MutationCount=trinucleotide()	
Bra2File0, Dec2Anc0, Anc2Dec0, Bra2MutC0, Go, TotalMutC, EdgeLs0 = Param.parse_config_file(Fol)
if Go=='y':
   print('branch and mutation count',Bra2MutC0,TotalMutC)
   Bra2MutCOri=Bra2MutC0   
   Cont='y'
   while Cont=='y':
        NewBra2OldBraID,NewBra2MutC = Format_MutationCount.MakeCombo(Dec2Anc0, Anc2Dec0, Bra2MutC0, EdgeLs0,MinMutC)
        if NewBra2MutC=={}: Cont='n'		
        else: Dec2Anc0, Anc2Dec0, Bra2MutC0, EdgeLs0=Format_MutationCount.UpBranchID(Dec2Anc0, Anc2Dec0, Bra2MutC0, EdgeLs0, NewBra2OldBraID,NewBra2MutC)
        print(NewBra2MutC)		
else: 
    summary_out=Go		
    print (summary_out)
if EdgeLs0==[]: 
    Go='n'
    summary_out='After combining short branches, there is no tree structure. So do not try to make an inference.'
    print(summary_out)
if Go=='y':
    print('combine small branches and make mutation count table')
    FinalBra2OldBraID={}
    FinalBra2MutC={}
    C=1
    for NBra in Bra2MutC0:
      if (NBra in Bra2MutCOri)==True: FBra=NBra
      else:	  
        FBra='C'+str(C)	
        C+=1
        FinalBra2OldBraID[FBra]=[NBra]			
      FinalBra2MutC[FBra]=Bra2MutC0[NBra]   
   	  
    Dec2Anc, Anc2Dec, Bra2SNVnum, EdgeLs=Format_MutationCount.UpBranchID(Dec2Anc0, Anc2Dec0, Bra2MutC0, EdgeLs0, FinalBra2OldBraID,FinalBra2MutC)
  	 
    Root,TipLs=Format_MutationCount.get_root_tip(Dec2Anc, Anc2Dec)
    Bra2AncSibDec=Format_MutationCount.get_ancsibdec(Dec2Anc, Anc2Dec)
    print('make input for PhyloSignare')
    Bra2Neigh={}
    Bra2File={}	
    out='New branch ID\tOriginal branch ID\tmutation count file\tNumber of mutations\n'
    for Bra in Bra2SNVnum:
        ComboCountCsv=ResDir+Bra+'_MutCount.csv'	
        if (Bra in FinalBra2OldBraID)==True: 
             Format_MutationCount.MakeMutCountFile(FinalBra2OldBraID[Bra][0], Bra2File0,ComboCountCsv)
             out+=Bra+'\t'+FinalBra2OldBraID[Bra][0]+'\t'+ComboCountCsv+'\t'+str(Bra2SNVnum[Bra])+'\n'			 
        else: 	
            shutil.copy2(Bra2File0[Bra], ComboCountCsv)					
            out+=Bra+'\t'+Bra+'\t'+ComboCountCsv+'\t'+str(Bra2SNVnum[Bra])+'\n'
        Bra2File[Bra]=ComboCountCsv		
        Neigh, TotMutC = Format_MutationCount.GetNeighboringBranch(Bra,Bra2SNVnum, Dec2Anc, Anc2Dec)	
        if TotMutC<100:
             NeighNall=[]	
             print(Bra, Neigh)			 
             for NbraLs in Neigh:
                if len(NbraLs)==1:			
                   NeighN, AA = Format_MutationCount.GetNeighboringBranch(NbraLs[0],Bra2SNVnum, Dec2Anc, Anc2Dec)				   
                   for Can in NeighN:
                      if Can!=[Bra]:
                           Can.append(NbraLs[0])					  
                           Can=list(set(Can))	
                           CanIn=[]
                           for CI in CanIn:
                              if CI!=Bra: CanIn.append(CI)						   
                           if CanIn!=[]:NeighNall.append(CanIn)						   
             Neigh+=NeighNall
        Bra2Neigh[Bra]=Neigh
    out+='#Tree\n'
    for E in EdgeLs:  out+=E+'\n'
    os.chdir(cwd)		
    SaveFile.GetOut(ResDir+'Summary.txt',out)	

ComboIDLs=[]
Unknown=[]
if Go=='y':	
    Bra2HitSig1={}
    for Bra in Bra2Neigh:
      print('compute signature activity 1, ',Bra)    	
      InF=Bra2File[Bra]
      if Bra2SNVnum[Bra]<5: pass	
      else: 
           SigTest=significanse_test(InF,Expected_SigLs,COS,TotSig,Rpath,CutSig)		   
           if os.path.exists(InF[:-4]+'_SSE.txt')!=True:		   		  
            Expo_All=SigTest.test_SSM(EstimateExposureRtemplate,InF[:-4]+'_SSE.txt',Method)	
            SaveFile.Clean(InF[:-4]+'_Rm*.txt')	
            shutil.copyfile(InF[:-4]+'_'+Method+'.txt', InF[:-4]+'_'+Method+'1.txt') 
            os.remove(InF[:-4]+'_'+Method+'.txt')  			
           RefitSig,OtherSigLs=SigTest.GetSignifSig1(InF[:-4]+'_SSE.txt', SigDif,BigSig)	 			   		   
           Bra2HitSig1[Bra]=RefitSig	
    BraCom2HitSig1={}
    BraCom2SSESig1={}	
    for Bra in Bra2Neigh: 
        InF=Bra2File[Bra]
        print('combine neighbors')		
        NeighLsLs=Bra2Neigh[Bra]	
        for NeighLs in NeighLsLs:
            print('make combo',Bra)
            CombFileLs=[InF]
            BraLs=[Bra]
            SNVnum=Bra2SNVnum[Bra]	         			
            for Neigh in NeighLs:
                CombFileLs.append(Bra2File[Neigh])
                BraLs.append(Neigh)
                SNVnum+=Bra2SNVnum[Neigh]				
            BraLs.sort()
            OutFileName=ResDir+''.join(BraLs)+'.csv'
            ComboIDLs.append(OutFileName)			
            print('combo-outfile name',OutFileName)
            print('pool input file list',CombFileLs)					            			
            Format_MutationCount.pool(CombFileLs,OutFileName)
            SigTest=significanse_test(OutFileName,Expected_SigLs,COS,TotSig,Rpath,CutSig)			
            if os.path.exists(OutFileName[:-4]+'_SSE.txt')!=True:			  
             Expo_All=SigTest.test_SSM(EstimateExposureRtemplate,OutFileName[:-4]+'_SSE.txt',Method)	
             SaveFile.Clean(OutFileName[:-4]+'_Rm*.txt')	
            RefitSig,OtherSigLs=SigTest.GetSignifSig1(OutFileName[:-4]+'_SSE.txt', SigDif,BigSig)	 	
            BraCom2SSESig1[Bra+'-'+'-'.join(map(str,NeighLs))]=RefitSig
    Bra2HitSig2={}
    Bra2SSESig2={}	
    for Bra in Bra2Neigh:
       if (Bra in Bra2HitSig1)!=True:	Bra2SSESig2[Bra]=[]	
       else:	
            print('collect signarures hit only when it is combo')	    	
            InF=Bra2File[Bra]
            SelSig=Bra2HitSig1[Bra]
            OriSig2SSEdif=SaveFile.GetSSE(InF[:-4]+'_SSE.txt')			
            NeighLsLs=Bra2Neigh[Bra]			
            CanLs=[]
            SSELs=[]
            SSE2SigLs={}			
            for NeighLs in NeighLsLs:
                 BNeighLs=NeighLs+[Bra]
                 BNeighLs.sort()
                 CSig2SSEdif=SaveFile.GetSSE(ResDir+''.join(map(str,BNeighLs))+'_SSE.txt')	         				 
                 CombSigLs=BraCom2SSESig1[Bra+'-'+'-'.join(map(str,NeighLs))]	
                 for CSig in CombSigLs:
                      if SelSig.count(CSig)==0:
                           if CanLs.count(CSig)==0: CanLs.append(CSig)
                           CSig2SSEdif=SaveFile.GetSSE(ResDir+''.join(map(str,BNeighLs))+'_SSE.txt')
                           CSSE=float(CSig2SSEdif['S'+str(CSig)])
                           if (CSSE in SSE2SigLs)!=True:
                                  SSE2SigLs[CSSE]=[]
                                  SSELs.append(CSSE)
                           SSE2SigLs[CSSE].append(CSig)								  		
            if CanLs==[] or (len(CanLs)==1 and SelSig==[]):			
                Bra2HitSig2[Bra]=SelSig
                Bra2SSESig2[Bra]=SelSig
                shutil.copyfile(InF[:-4]+'_SSE.txt', InF[:-4]+'_SSE1.txt')
                shutil.copyfile(InF[:-4]+'_'+Method+'1.txt', InF[:-4]+'_'+Method+'2.txt')				 
            else:
                SSELs.sort(reverse=True)
                AllLs=[]
                Pos=0				
                while len(AllLs)<10 and Pos<len(SSELs):
                    SSE=SSELs[Pos]
                    New=SSE2SigLs[SSE]
                    New=list(set(New))					
                    AllLs+=New
                    AllLs=list(set(AllLs))					
                    Pos+=1									
                AllLs+=SelSig
                AllLs=list(set(AllLs))			   
                AllLs.sort()		   
                AllLs_str=','.join(map(str, AllLs))	
                SigTest=significanse_test(InF,AllLs_str,COS,TotSig,Rpath,CutSig)				   
                if os.path.exists(InF[:-4]+'_SSE1.txt')	!=True:				
                 Expo_All=SigTest.test_SSM_limit1(EstimateExposureRtemplate,InF[:-4]+'_SSE1.txt',Method,AllLs)					
                 SaveFile.Clean(InF[:-4]+'_Rm*.txt')	
                 shutil.copyfile(InF[:-4]+'__'+Method+'.txt', InF[:-4]+'_'+Method+'2.txt') 
                 os.remove(InF[:-4]+'__'+Method+'.txt')  					 
                GoodCanSig,OtherSigLs=SigTest.GetSignifSig1(InF[:-4]+'_SSE1.txt', 0.05,BigSig)	 
                Bra2HitSig2[Bra]=SelSig+GoodCanSig+OtherSigLs					
                Bra2SSESig2[Bra]=SelSig+GoodCanSig
    print('make output file')		   
    for Bra in Bra2Neigh: 
           if os.path.exists(ResDir+Bra+'_MutCount_'+Method+'2.txt')==True and Bra2SNVnum[Bra]>=OutPutCut:		  	  
            InF=Bra2File[Bra]
            Supp=Bra2SSESig2[Bra]	
            Sig2SSEdif=SaveFile.GetSSE(InF[:-4]+'_SSE1.txt')	
            AA,Sig2Est=SaveFile.GetHitSig(ResDir+Bra+'_MutCount_'+Method+'2.txt',0.01)
           else: summary_out+='PhyloSigFinder did not make an inference because the number of mutations was small: '+Bra+'\n'				
    print('test signature change between branches')
    Try=1
    Bra2SSESig3={}	
    Conti='y'	
    if Conti=='y':
          for Bra in Bra2Neigh:
           if Bra2SSESig2[Bra]==[] and Bra2SNVnum[Bra]<OutPutCut: Bra2SSESig3[Bra]=[]
           else: 		   
            print('neighbor is SSM signature + >0.01 activity in combo', Bra,Bra2SSESig2[Bra])	
            NeighLsLs=Bra2Neigh[Bra]	
            OriSSESig=Bra2SSESig2[Bra]	
            InF=Bra2File[Bra]				
            OriHitSigLs,OriSig2Est=SaveFile.GetHitSig(ResDir+Bra+'_MutCount_'+Method+'2.txt',0.05)	
            OriSig2SSEdif=SaveFile.GetSSE(InF[:-4]+'_SSE1.txt')			
            OriSig2SSEdif0=SaveFile.GetSSE(InF[:-4]+'_SSE.txt')			
            NeiSSESigLs=[]
            ComboHit=[]		
            GoodLs=[]			
            for NeighLs in NeighLsLs: 
                 NBra=NeighLs[0]
                 BNeighLs=NeighLs+[Bra]
                 BNeighLs.sort()
                 CHitSigLs=SaveFile.GetSSE(ResDir+''.join(map(str,BNeighLs))+'_SSE.txt')				 
                 ComboHit+=CHitSigLs
                 for NBra in NeighLs:
                   if (NBra in Bra2HitSig1)==True:				 
                       NeiSSESig=Bra2HitSig1[NBra]
                       NeiSSESigLs+=NeiSSESig
            NeiSSESigLs=list(set(NeiSSESigLs))
            ComboHit=list(set(ComboHit))
            CanLs=[]
            for NS in NeiSSESigLs:
                 if ComboHit.count('S'+str(NS))!=0 and OriSSESig.count(NS)==0 and OriHitSigLs.count(NS)!=0 :				 
                        GoodLs.append(NS)				 				 
                 elif ComboHit.count('S'+str(NS))!=0 and OriSSESig.count(NS)==0: CanLs.append(NS)
                 elif OriSSESig==[]: CanLs.append(NS)
            if  OriSSESig==[] and CanLs==[]: 
                for i in ComboHit: 
                   if i!='SAll':CanLs.append(int(i.replace('S','')))              				 	   
            if CanLs==[]:
                 Bra2SSESig3[Bra]=OriSSESig+GoodLs
            else:				 
             CanLs+=OriSSESig+GoodLs	
             CanLs=list(set(CanLs))			 
             CanLs.sort()				 
             Cand_SigLs=','.join(map(str,CanLs))			
             Estimate_Exposure=Exposure_computer(InF,Cand_SigLs,COS,TotSig,EstimateExposureRtemplate,Rpath,CutSig)	
             Expo, SSE,SSE_per_MutType=Estimate_Exposure.estimate_exposureR(InF[:-4]+'_'+Method+'3.txt')
             RefitSig=Estimate_Exposure.Get_largeSig_all(Expo,0.05) 			 
             New=RefitSig+OriSSESig+GoodLs
             New=list(set(New))					 
             Bra2SSESig3[Bra]=New
        #  print (Bra2SSESig3)	  
          Bra2SSESig4={}
          print('sequential loss and gain')	
          Bra2CorTestSig={}	
          for Bra in Bra2Neigh:
           if Bra not in Bra2CorTestSig: 			  
               Bra2CorTestSig[Bra]={}		  
           if Bra2SSESig3[Bra]==[]: 	   
             Bra2CorTestSig[Bra]['Rec-Evo']=[]
             if 'Amb-Evo' not in Bra2CorTestSig[Bra]:			 
                 Bra2CorTestSig[Bra]['Amb-Evo']=[]
           else: 		  
            RecLs=[]
            AmbLs=[]
            CanLs=[]			
            OriSSESig=Bra2SSESig3[Bra]			
            if Bra!=Root and TipLs.count(Bra)==0:
                InF=Bra2File[Bra]			
                AncSibDec=Bra2AncSibDec[Bra]				
                AncSSESig=Bra2SSESig3[AncSibDec['Anc']]								
                SibSSESig=[]
                SibLs=AncSibDec['Sib']
                for Sib in SibLs:
                  	SibSSESig+=Bra2SSESig3[Sib]	
                DecLs=AncSibDec['Dec']
                DecSSESig=[]				
                for Dec in DecLs:
                  	DecSSESig+=Bra2SSESig3[Dec]  
                DecSSESigSet=list(set(DecSSESig))			
                for Dsig in DecSSESigSet:
                  if OriSSESig.count(Dsig)==0:				
                    if SibSSESig.count(Dsig)!=0 or AncSSESig.count(Dsig)!=0:
                       MissC=len(DecLs)-DecSSESigSet.count(Dsig)+len(SibLs)-SibSSESig.count(Dsig)+1-AncSSESig.count(Dsig)
                      # print (Bra,Dsig,MissC,DecLs,DecSSESigSet,SibLs,SibSSESig,AncSSESig)
                       if MissC<=1:					                       					   
                         CanLs.append(Dsig)		
                       else:
                         if AncSSESig.count(Dsig)!=0: 
                            if AncSibDec['Anc'] not in Bra2CorTestSig:Bra2CorTestSig[AncSibDec['Anc']]={'Amb-Evo':[]}						
                            Bra2CorTestSig[AncSibDec['Anc']]['Amb-Evo'].append(Dsig)
                         if DecSSESig.count(Dsig)!=len(DecLs): 
                            for Dec in DecLs:
                                if Bra2SSESig3[Dec].count(Dsig)!=0:						 
                                   if Dec not in Bra2CorTestSig:Bra2CorTestSig[Dec]={'Amb-Evo':[]}			
                                   Bra2CorTestSig[Dec]['Amb-Evo'].append(Dsig)
                         if SibSSESig.count(Dsig)!=len(SibLs): 
                            for Sib in SibLs:
                                 if Bra2SSESig3[Sib].count(Dsig)!=0:  						 
                                   if Sib not in Bra2CorTestSig:Bra2CorTestSig[Sib]={'Amb-Evo':[]}							
                                   Bra2CorTestSig[Sib]['Amb-Evo'].append(Dsig)								   
            if len(CanLs)==0:
               Bra2SSESig4[Bra]=OriSSESig			
            else: 
              AllLs=CanLs+OriSSESig
              AllLs=list(set(AllLs))			  
              if len(AllLs)<=1: Bra2SSESig4[Bra]=OriSSESig
              else:			  
                AllLs.sort()				 
                Cand_SigLs=','.join(map(str,AllLs))		
                Estimate_Exposure=Exposure_computer(InF,Cand_SigLs,COS,TotSig,EstimateExposureRtemplate,Rpath,CutSig)	
                Expo, SSE,SSE_per_MutType=Estimate_Exposure.estimate_exposureR(InF[:-4]+'_'+Method+'4.txt')
                RefitSig=Estimate_Exposure.Get_largeSig_all(Expo,CutSig) 
                for CS in CanLs:
                    if RefitSig.count(CS)!=0: RecLs.append(CS)
                    else: AmbLs.append(CS)					
                New=OriSSESig+RefitSig
                New=list(set(New))				
                Bra2SSESig4[Bra]=New	
            Bra2CorTestSig[Bra]['Rec-Evo']=RecLs
            if 'Amb-Evo' not in Bra2CorTestSig[Bra]:Bra2CorTestSig[Bra]['Amb-Evo']=[]				
            Bra2CorTestSig[Bra]['Amb-Evo']+=AmbLs
         # print (Bra2CorTestSig)	 
          print('make output file')			   
          AllResIn='Signature\tRelative activity ('+Method+')\tACS\tchi test\tBranch\tflag\n'	
          m2pval=SaveFile.readChiTa(ChiTa)		  
          for Bra in Bra2Neigh: 
           if os.path.exists(ResDir+Bra+'_MutCount_'+Method+'2.txt')==True and Bra2SNVnum[Bra]>=OutPutCut:		    
            InF=Bra2File[Bra]
            if (Bra in Bra2SSESig4)==True:			
             Supp=Bra2SSESig4[Bra]	
             Supp=list(set(Supp))			 
             Sig2SSEdif=SaveFile.GetSSE(InF[:-4]+'_SSE1.txt')	
             AA,Sig2Est=SaveFile.GetHitSig(ResDir+Bra+'_MutCount_'+Method+'2.txt',0.01)						
             Amb=Bra2CorTestSig[Bra]['Amb-Evo']			  
             SaveFile.PhyloSigFinder_output3(Sig2SSEdif,Supp,Amb,Sig2Est,TotSig,ResDir+Bra+'_MutCount_PhyloSigFinder.txt') 		 
             print ('final SSE test and activity inference',Bra)			 
             SigTest=significanse_test(Bra2File[Bra],','.join(map(str, Supp)),COS,TotSig,Rpath,CutSig)			
             Expo_All=SigTest.test_SSM(EstimateExposureRtemplate,Bra2File[Bra][:-4]+'final_SSE.txt',Method)	
             Sig2SSEdifFinal=SaveFile.GetSSE(Bra2File[Bra][:-4]+'final_SSE.txt')
             Expo_All0=SaveFile.list_exposure(Bra2File[Bra][:-4]+'final_'+Method+'_All.txt')			 
             Chi_all,df=SigTest.compute_chi(Expo_All0)	
             m=(df-1)/2	
           #  print ('m',m,'chi-all',Chi_all)			 
             for S in Supp:
                # print ('rm sig',S)
                 if Amb.count(S)!=0: Flag='?'
                 else: Flag=''					 
                 EstFile=Bra2File[Bra][:-4]+'_Rm'+str(S)+'.txt'		
                 if os.path.exists(EstFile)!=True: AllResIn+='S'+str(S)+'\tActive\tNA\tNA\t'+Bra+'\t'+Flag+'\n'
                 else: 				 
                  RmExpLs=SaveFile.list_exposure(EstFile) 
                  Chi_Rm,AA=SigTest.compute_chi(RmExpLs)
                  if Chi_Rm<=Chi_all: P='>0.05'
                  else: 
                      DeltaChi=(Chi_Rm-Chi_all)/2
                    #  print ('DeltaChi (rm-all)',DeltaChi,Chi_Rm,Chi_all)					  
                      P=SigTest.chitest(DeltaChi,m,m2pval)			  
                  AllResIn+='S'+str(S)+'\t'+str(Expo_All[S-1])+'\t'+str(Sig2SSEdifFinal['S'+str(S)])+'\t'+P+'\t'+Bra+'\t'+Flag+'\n'					                       					  
                  SaveFile.Clean(EstFile)	
            # print (AllResIn)			 
             SaveFile.GetOut(ResDir+'PhyloSignare.txt',AllResIn)			
   		
for Combcsv in ComboIDLs: 
    SaveFile.Clean(Combcsv)
    SaveFile.Clean(Combcsv[:-4]+'_*.txt')
ParamSig.map_signature_ID(ResDir+'PhyloSignare.txt')
SaveFile.GetOut(ResDir+'Note.txt',summary_out)

if FigOut=='y':
    print ('Create GraphViz')
    ColorFile=sys.argv[3]
    Color=open(ColorFile,'r').readlines()
    Sig2Color={}
    for i in Color:
        i=i.strip().split('\t')
        Sig2Color[i[0]]=i[1]
   # from output.Node import Node
    from graphviz import Source
    from pathlib import Path
    os.chdir(ResDir)		
    class Node:
        tree_age = 1
    
        def __init__(self, name, number_of_mutations):
            self.parent = None
            self.age = -1
            self.isNumbered = False
            self.children = []
            self.name = name
            self.number_of_mutations = number_of_mutations
            self.isChild = False
            self.numbered_name = "Normal"
    
        def __str__(self):
            return ', label = "' + 'B' + str(self.age) + ': ' + self.number_of_mutations + '"'
    
        # Looks for phylosigfinder for appropriate node and determines how many presences there are
        @staticmethod
        def read(node, path):
    
            # Opens phylosigfinder file
    
            if Path(path / Path(node.name + '_MutCount_PhyloSigFinder.txt')).exists():    
                phylosigfinder = open(path / Path(node.name + '_MutCount_PhyloSigFinder.txt'))
                phylosigfinder = phylosigfinder.readlines()
                presences = []
            
                        # Remove the first line
                phylosigfinder.pop(0)
            
                        # Split each line into an array
                for i in range(len(phylosigfinder)):
                    phylosigfinder[i] = phylosigfinder[i].split()
            
                        # Detects which signatures are present
                for line in phylosigfinder:
                    if line[1] == '1':
                        presences.append(line[0])
                return presences
            else:
                return []
    
        def add_child(self, node):
            node.parent = self
            self.children.append(node)
    
        def set_age(self):
            self.age = Node.tree_age
            Node.tree_age += 1
            self.isNumbered = True
            print(self.name + " age: " + str(self.age))
    
 	
    def update(tree, nodes):
        node_names = []
        new_tree = []
    
        # Creates a list of node names, following same order
        for node in nodes:
            node_names.append(node.name)
    
        # Reads tree to add children to node
        for branch in tree:
    
            # Assigns child node to parent
            parent = nodes[node_names.index(branch[:branch.index('-')])]
            child = nodes[node_names.index(branch[branch.index('>') + 1:])]
            child.isChild = True
    
            parent.add_child(child)
    
            #Prepares new tree
            new_tree.append([parent, child])
    
        # Assign ages
        current_node = None
        for node in nodes:
            # Find parent node
            if node.isChild:
                continue
            # Assign parent node
            current_node = node
    
            # Assign parent node age to 1
            node.set_age()
    
        while True:
            # Iterate over tree
            while len(current_node.children) > 0:
                print("Current node: " + current_node.name)
                low_index = 0
    
                # Find the child with least children
                for i in range(len(current_node.children)):
                    print(current_node.children[i].name + ":" + str(len(current_node.children[i].children)) + "-" + str(current_node.children[i].isNumbered))
                    if len(current_node.children[i].children) <= len(current_node.children[low_index].children) or not current_node.children[i].isNumbered:
                        if i == low_index:
                            print("comparing same node, moving on")
                        else:
                            print(current_node.children[i].name + " has less children than " + current_node.children[low_index].name)
                            if current_node.children[i].isNumbered:
                                print(current_node.children[i].name + " is already numbered, moving on")
                            else:
                                print(current_node.children[i].name + " is not numbered yet")
                                low_index = i
    
                print("Lowest child: " + current_node.children[low_index].name)
                # Set age and move down the tree
                if current_node.children[low_index].isNumbered:
                    print("Child already numbered")
                    break
                current_node.children[low_index].set_age()
                current_node = current_node.children[low_index]
    
            # Back up tree
            print("End of line, setting current node to: " + current_node.parent.name)
            current_node = current_node.parent

            all_numbered = True
            for node in nodes:
                if not node.isNumbered:
                    print(node.name + " is not numbered, continuing")
                    all_numbered = False
            if all_numbered:
                break
        print("all nodes numbered")
        Node.tree_age = 1
    
        # Rename each node based on age
        for node in nodes:
            node.numbered_name = "N" + str(node.age)
    
        return new_tree
    

    def write(path):
        # Read the Summary.text and create an array containing each line
        summary = open('Summary.txt')
        summary = summary.readlines()
        nodes = []
        tree = []
    
        # Remove the first line
        summary.pop(0)
    
        # Find #Tree and split each line into an array
        temp = -1
        for i in range(len(summary)):
            if summary[i] == '#Tree\n':
                temp = i
            summary[i] = summary[i].split()
    
        # Modify summary and tree
        tree = summary[temp + 1:]
        for i in range(len(tree)):
            tree[i] = str(tree[i]).strip("[']")
        summary = summary[0:temp]
    
    
    
        # Fill out nodes
        for i in range(len(summary)):
            nodes.append(Node(summary[i][0], summary[i][3]))
    
        # Creating the Output file
        output_file = open("output.gv", "a")
    
        # Writing the opening to the output file
        output_file.write('digraph ' + '"' + str(sys.argv[1]) + '"' + '{\n\nrankdir="LR";\n')
        output_file.close()
    
        # For each branch-pair in tree, add it to output
        output_file = open("output.gv", "a")
    
        # Updates node and tree
        tree = update(tree, nodes)
    
        # Adds the normal pair to the tree
        for node in nodes:
            if node.numbered_name == 'N1':
                tree.append([Node('Normal', -1), node])
    
        # Reads the last 2 letters of each element in tree and reads that node
        for pair in tree:
            # Reads the node given
            presences = Node.read(pair[1], path)
            color = ""
    
            # Fills out draw_line code
    
            # Case file does not exist
            print(presences)
            if len(presences) == 0:
                color = 'black'
                label = str(pair[1])
                draw_line = ' [dir = none, penwidth = 2, color = ' + color + label + ']\n'
                output_file.write(pair[0].numbered_name + '->' + pair[1].numbered_name + draw_line)
                
            else:
    
                colorlist=[v for k, v in Sig2Color.items() if k in presences]
                colorname = ':white:'.join(colorlist)
                label = str(pair[1])
    
                draw_line = ' [dir = none, penwidth = 2, color = ' + '"' + colorname + '"' + label + ']\n'
                output_file.write(pair[0].numbered_name + '->' + pair[1].numbered_name + draw_line)
    
        output_file.close()
    
        # Writing the closing to the output file
        output_file = open("output.gv", "a")
        output_file.write('labelloc="t";\nlabel="' + str(sys.argv[1]) + '"; \n\n}\n')
        output_file.close()
    
        print("Added " + path)


    open(Path.cwd() / Path('output.gv'), 'w').close()
    currentpath = os.path.dirname(os.path.abspath(__file__))
    print(currentpath)

    write(currentpath)
    outputgv = os.path.join(currentpath, 'output.gv')
    s = Source.from_file(outputgv)
    s.view()		
os.chdir(cwd)			
for Bra in Bra2File:
    FileLs=glob.glob(Bra2File[Bra][:-4]+'*.txt')
    for File in FileLs:
            SaveFile.Clean(File)	
			
