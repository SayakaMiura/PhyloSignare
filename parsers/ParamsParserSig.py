
import os

class ParamsParserSig(object):
    
    def __init__(self,ContFile):        
        self.ContFile=ContFile
        self.Method2TMP={'QP':'QP-run.r'}
       # self.Method2TMP['SA']='SA-run.r'
        self.Method2TMP['MutCon']='MutCon-run.r'
        self.Method2TMP['MutPat']='MutPat-run.r'
        #self.Method2TMP['sigfit']='sigfit-run.r'
        #self.Method2TMP['sigL']='C:\\Users\\tuf78332\\Desktop\\mutational-signature\\PhyloSig\\PhyloSig1\\sigL-run.r'
        self.Method2TMP['dSig']='dSig-run.r'        
       
    def get_parameter_setting(self):
         Con=open(self.ContFile,'r').readlines()
      #   print Con		 
         COS='Signature table is not given'
         Rpath='Rpath is not given'
         self.Can2SigLs={}#'signature list file is not given'		 
       #  CutSig='cut off value for signature is not given'
         self.Method='base method for signature refitting is not given'
         self.SumOut=''	
         ID=''		 
         SigLsIn=[]        
         for i0 in Con:
            i=i0.strip().split('\t')
           # print i[0],':::',i[1]			
            if i[0]=='BaseMethod': self.Method=i[1]
            elif i[0]=='Signature': 
              #  print i[1].strip()			
                if int(i[1])==2: 
                    COS='cosmicSig.csv'
                    Translate='cosmicSigv2IDlist.txt'					
                elif int(i[1])==3: 
                    COS='cosmicSigv3.csv'	
                    Translate='cosmicSigv3IDlist.txt'						
                else: 
                   COS='it should be 2 or 3'	
                   self.SumOut+=COS	
               # print 'h',COS,int(i[1])				   
            elif i[0]=='Signature List':
               Translate=open(Translate,'r').readlines()
               self.In2Out={}
               self.Out2In={}				   
               for Ti in Translate:
                       Ti=Ti.strip().split('\t')
                       self.In2Out[Ti[1]]=int(Ti[0].replace('S',''))
                       self.Out2In[Ti[0]]=Ti[1]	 			
               if i[1].strip()=='all' or i[1].strip()=='All':
                    if COS=='cosmicSig.csv': self.Can2SigLs['All']=list(range(1,31))
                    elif COS=='cosmicSigv3.csv': self.Can2SigLs['All']=list(range(1,66))
                   # print COS					
                    SigLsIn=self.Can2SigLs['All']	
                    ID='All'					
               else: 
				   
                   TargetLs=i[1].split(',')
                 #  print TargetLs,In2Out,i				   
                   for Target in TargetLs:
                       if (Target in self.In2Out)!=True: self.SumOut+='Signature ID is not found in COSMIC. Please correct: '+Target+'\n'
                       else: SigLsIn.append(self.In2Out[Target])
                       SigLsIn=list(set(SigLsIn))					   
                  # SigLsIn=[int(i) for i in i[1].split(',')]
                  # self.Can2SigLs['Part']=SigLsIn				   
            elif i[0]=='Signature ID':
               if ID=='': self.Can2SigLs[i[1].strip()]=SigLsIn
            elif i[0]=='Rpath': Rpath='\"'+i[1]+'\"'
            else:
               self.SumOut+='incorrect value in the Control.txt'+': '+i0	+'\n'
               print('incorrect value in the Control.txt'+': '+i0)	
         self.SumOut+='parameter setting\nSignature table: '+COS+'\n'+'Rpath: '+Rpath+'\n'+'Base method: '+self.Method+'\n'		
        # print self.Out2In
       #  open('A','r').readlines()		 
         return COS,Rpath,self.Method,self.Can2SigLs,SigLsIn,self.Method2TMP[self.Method]
    def map_signature_ID(self,File):
        OutFile=File[:-4]+'-clean.txt'		
        File1=open(File,'r').readlines()
        out=File1[0].replace('\tACS\t','\tiS\t')
        File1=File1[1:]
     #   print 'h',self.Out2In,File		
        for i in File1:
           Sid=i.split('\t')[0].split(' ')[0]
         #  print Sid		   
           if Sid=='SAll': out+=i
           else:		   
             In=self.Out2In[Sid]+i[len(Sid):]		   
             out+= In	
        OutF=open(File,'w')
        OutF.write(out)
        OutF.close()		
	
    def summary_signature_input_list(self,File):
           ID2COS={}
           COS2ID={}
           SigLs=[]
           File=open(File,'r').readlines()
           c=0		   
           for i in File:
               c+=1				   
               i=i.strip().split('\t')
               ID=int(i[0].replace('S',''))
               if ID!=c:
                   print('ID is incorrect: '+File+str(ID))
                   self.SumOut+='ID is incorrect: '+File+str(ID)
              				   
               COS=i[1]
               PreAbs=int(i[2])
               ID2COS[ID]=COS
               COS2ID[COS]=ID
               if PreAbs==1: SigLs.append(ID)
	   
           return ID2COS,COS2ID,SigLs,c			   
    def get_signature_id(self):
         ID2COS={}
         COS2ID={}
         Expected_SigLsDic={}		 
         if self.Can2SigLs=='signature list file is not given': 
             self.SumOut+='signature list file is not given\n'
             return {},{},{}			 
         else:
             for ID in self.Can2SigLs:
                 File=self.Can2SigLs[ID]
                 if os.path.exists(File)!=True: 
                    self.SumOut+='Signature list file does not exist: '+File+'\n'
                    print('Signature list file does not exist: '+File+'\n')
                 else:
                    ID2COS,COS2ID,SigLs,TotSig=self.summary_signature_input_list(File)	
                    Expected_SigLsDic[ID]=	SigLs				
			 
         		 
         return ID2COS,COS2ID,Expected_SigLsDic,TotSig	
    def get_summary_out(self):
         return self.SumOut	
    def get_Rtemplate(self):
         if (self.Method in self.Method2TMP)!=True: 
              print('base methods should be QP,SA,dSig,MutPat,or MutCon: '+self.Method)
              self.SumOut+=	'base methods should be QP,SA,dSig,MutPat,or MutCon: '+self.Method+'\n'		  
         return self.Method2TMP[self.Method]	