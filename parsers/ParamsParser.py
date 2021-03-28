
import os

class ParamsParser(object):
    
    def __init__(self,Input_file):        
        self.InFile=Input_file
        
       
    def CountMutation(self,csv):
         csv=open(csv,'r').readlines()
         Mcount=0
         for i in csv:
            Mcount+=int(i.split(',')[1])
         return Mcount	   
    def parse_config_file(self,RootDir):
      if os.path.exists(self.InFile)!=True:
         print('input file does not exists: ', self.InFile)
         return '', '', '', '', 'n' ,'',''		 
      else:		 
         File=open(self.InFile,'r').readlines()[1:]
         Tree='n'
         Go='y'	
         Dec2Anc={}
         Anc2Dec={}
         Branch2McouFile={}
         Branch2Mcou={}
         TreeLine=[]	
         Tot=0	
         for i in File:
             if Tree=='y':
                 TreeLine.append(i.strip())		
                 i=i.strip().split('->')
                 Anc=i[0]
                 Dec=i[1]
                 if (Anc in Anc2Dec)!=True: Anc2Dec[Anc]=[]
                 Anc2Dec[Anc].append(Dec)
                 Dec2Anc[Dec]=Anc
             elif i[0]=='#': Tree='y'
             else:
                i=i.strip().split('\t')
                Bra=i[0]
                McouFile=RootDir+i[1]#.replace('\\','\\\\')		   
                if os.path.exists(McouFile)!=True:
                    print('no file',McouFile)
                    Go='n'
                else:			   
                    Mc=self.CountMutation(McouFile)		
                    Branch2Mcou[Bra]=Mc
                    Branch2McouFile[Bra]=McouFile
                    Tot+=Mc	
         if Tot<100: 
            print('the number of total mutations is <100, so do not try inference',Tot)
            Go='the number of total mutations is <100, so do not try inference: '+str(Tot)	   
         return Branch2McouFile, Dec2Anc, Anc2Dec, Branch2Mcou, Go ,Tot,TreeLine	