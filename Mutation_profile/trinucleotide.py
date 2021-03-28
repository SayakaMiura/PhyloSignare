from output.file_writer import file_writer

class trinucleotide():
     def Rename(self,Dec2Anc0, Anc2Dec0, Bra2MutC0, EdgeLs0, NewBra2OldBraID):
        Dec2Anc1={}
        Anc2Dec1={}
        Bra2MutC1={}
        EdgeLs1={}		
     def GetTip(self,Dec2Anc,Anc2Dec):
        Tips=[]
        for Dec in Dec2Anc:
           if (Dec in Anc2Dec)!=True: Tips.append(Dec)
        return Tips	
     def get_root_tip(self,Dec2Anc, Anc2Dec):
         self.Root=''
         self.TipLs=[]
         print(Anc2Dec)		 
         for Dec in Dec2Anc:
           #print(Dec,Dec in Anc2Dec)
           if (Dec in Anc2Dec)!=True: self.TipLs.append(Dec)
         for Anc in Anc2Dec:
            if (Anc in Dec2Anc)!=True: self.Root=Anc		 
         print(self.TipLs)           				 
         return self.Root,self.TipLs	
     def get_ancsibdec(self,Dec2Anc, Anc2Dec):	
         Bra2Anno={}
         for Dec in Dec2Anc:
              Anc=Dec2Anc[Dec]
              SibLs0=Anc2Dec[Anc]
              SibLs=[]
              for Sib in SibLs0:
                 if Sib!=Dec:
                     SibLs.append(Sib)
              DecLs=[] 					 
              for Dec1 in Dec2Anc:
                 if Dec2Anc[Dec1]==Dec: DecLs.append(Dec1)
              Bra2Anno[Dec]={'Anc':Anc,'Sib':SibLs,'Dec':DecLs}
         return Bra2Anno			  
     def MakeAncCombo(self,IniBra,Dec2Anc, Anc2Dec,MutC,Bra2MutC,MinM): 
         print('make combo',IniBra,Dec2Anc,Anc2Dec)			 
         IniAnc=Dec2Anc[IniBra]
         Good='y'
         NewComb=[]		 
         if (IniBra in Anc2Dec)==True:		 
            if len(Anc2Dec[IniBra])>1: 
               Good='n'
               NewComb=[]
               NewMutCombC=MutC	
         if Good=='y':			   
          if len(Anc2Dec[IniAnc])>1: 
            Good='n'
            NewComb=[]
            NewMutCombC=MutC	   
          elif NewComb.count(IniAnc)==0:     
            Good='y'	
            NewComb=[IniAnc]
            NewMutCombC=MutC+Bra2MutC[IniAnc]
     
         while Good=='y' and NewMutCombC<MinM: #add up
            if (IniAnc in Dec2Anc)!=True: Good='n'
            else:
               IniAnc=Dec2Anc[IniAnc]
               Decs=Anc2Dec[IniAnc]
               if len(Decs)>1: Good='n'
               elif NewComb.count(IniAnc)==0:
                  NewComb.append(IniAnc)
                  NewMutCombC+=Bra2MutC[IniAnc]	
         Good1='y'
         while Good1=='y' and NewMutCombC<MinM: #add down
            if (IniBra in Anc2Dec)!=True: Good1='n'
            else:
               Decs=Anc2Dec[IniBra]
               if len(Decs)>1: Good1='n'
               elif NewComb.count(IniBra)==0:
                  IniBra=Decs[0]		  
                  NewComb.append(IniBra)
                  NewMutCombC+=Bra2MutC[IniBra]		
         print('ancestors added',NewComb,NewMutCombC)
         return NewComb,NewMutCombC	   
     def MakeCombo(self,Dec2AncOri, Anc2DecOri, Bra2MutCOri, EdgeLsOri,MinM): 
         BadBraLs=[]
         for B in Bra2MutCOri:
             if Bra2MutCOri[B]<MinM:BadBraLs.append(B)
         print('bad branch list',BadBraLs)		
         TipLs=self.GetTip(Dec2AncOri,Anc2DecOri)
         print('fix small tip branches, tips',TipLs)	
         NewBadBraLs=[]	
         Done=[]
         OldBra2NewBraID={}
         NewBra2OldBraID={}	
         NewBra2MutC={}	
         for Tip in TipLs: 
            TarTipMuc=Bra2MutCOri[Tip]	
            if BadBraLs.count(Tip)!=0 and Done.count(Tip)==0:
               print('number of mutations is small, so try to make combo',Tip,Bra2MutCOri[Tip]) 	   
               Anc=Dec2AncOri[Tip]
               Sibs=Anc2DecOri[Anc]
               TipSib2Mc={}
               Mc2TipSib={}			   
               for T in Sibs:
                  if TipLs.count(T)!=0 and T!=Tip:
                       Mc=Bra2MutCOri[T]			 
                       TipSib2Mc[T]=Mc
                       if (Mc in Mc2TipSib)!=True:Mc2TipSib[Mc]=[]
                       Mc2TipSib[Mc].append(T) 				  
               print('tip siblings',Mc2TipSib)	
               CombLs=[Tip]
               MutCombC=TarTipMuc		  
               if Mc2TipSib!={}: 
                   McLs=list(Mc2TipSib.keys())	
                   McLs.sort()
                   print(McLs)			  
     
                   LsTot=len(McLs)
                   LsP=0			  
                   while MutCombC<MinM and LsP<LsTot: 
                       SibLs=Mc2TipSib[McLs[LsP]]
                       for i in SibLs:
                           MutCombC+=Bra2MutCOri[i]
                           CombLs.append(i)
                       LsP+=1					  
                   Done+=CombLs
                   print('combo after adding siblings',CombLs,MutCombC)
               if MutCombC<MinM: ##Add Anc
                      AncComb,MutCombC=self.MakeAncCombo(Tip,Dec2AncOri, Anc2DecOri,MutCombC,Bra2MutCOri,MinM)	
                      Done+=AncComb
                      CombLs+=AncComb
               if len(CombLs)>1:
                  NewID=''.join(CombLs)		  
                  OldBra2NewBraID[Tip]=NewID
                  NewBra2OldBraID[NewID]=CombLs	
                  NewBra2MutC[NewID]=MutCombC				  
         print('fix small intermediate branches')			 
         for BadBra in BadBraLs: ##fix small intermediate branches
            if Done.count(BadBra)==0 and TipLs.count(BadBra)==0:
                  print('h',TipLs,BadBra,Anc2DecOri)			
                  Bottom=BadBra
                  More='y'
                  while More=='y':
                      Decs=Anc2DecOri[Bottom]
                      if len(Decs)>1: More='n'
                      else:
                         if Bra2MutCOri[Decs[0]]>=MinM: More='n'
                         else: 
                             Bottom=Decs[0]
                             if (Bottom in Anc2DecOri)!=True: More='n'							 
                  print('examine inter ',Bottom)
                  CombLs=[Bottom]			 
                  AncComb,MutCombC=self.MakeAncCombo(Bottom,Dec2AncOri, Anc2DecOri,Bra2MutCOri[Bottom],Bra2MutCOri,MinM)	
                  Done+=AncComb
                  CombLs+=AncComb
                  if len(CombLs)>1:
                      NewID=''.join(CombLs)		  
                      OldBra2NewBraID[Tip]=NewID
                      NewBra2OldBraID[NewID]=CombLs	
                      NewBra2MutC[NewID]=MutCombC 
         print(NewBra2OldBraID,'\n',NewBra2MutC)				 
         return NewBra2OldBraID,NewBra2MutC
     def UpBranchID(self,Dec2AncOri, Anc2DecOri, Bra2MutCOri, EdgeLsOri, NewBra2OldBra,NewBra2MutC):
         Dec2AncNew={}
         Anc2DecNew={}
         Bra2MutCNew={}
         EdgeLsNew=[]
         Old2New={}
     	
         for New in 	NewBra2OldBra:
             OldBraLs=NewBra2OldBra[New]
             for Old in OldBraLs:
                 Old2New[Old]=New
         for Dec in Dec2AncOri:
             Anc=Dec2AncOri[Dec]
             if (Anc in Old2New)==True: Anc=	Old2New[Anc]
             if (Dec in Old2New)==True: Dec=	Old2New[Dec]		
             if Dec!=Anc: Dec2AncNew[Dec]=Anc
         for Anc in Anc2DecOri:
             DecLs=Anc2DecOri[Anc]
             if (Anc in Old2New)==True: Anc=	Old2New[Anc]
             NewDecLs=[]		
             for Dec in DecLs:		
                 if (Dec in Old2New)==True: Dec=	Old2New[Dec]	
                 if Dec!=Anc: NewDecLs.append(Dec)
             NewDecLs=list(set(NewDecLs))			
             if NewDecLs!=[]: Anc2DecNew[Anc]=NewDecLs    		
         for Bra in Bra2MutCOri:
             MutC=Bra2MutCOri[Bra]	
             if (Bra in Old2New)==True: 
                   Bra=	Old2New[Bra]
                   MutC=	NewBra2MutC[Bra]		  
             Bra2MutCNew[Bra]=MutC
         for Edge in EdgeLsOri:
             Anc=Edge.split('->')[0]
             Dec=Edge.split('->')[1]
             if (Anc in Old2New)==True: Anc=	Old2New[Anc]
             if (Dec in Old2New)==True: Dec=	Old2New[Dec]
             if Anc!=Dec: EdgeLsNew.append(Anc+'->'+Dec)
         EdgeLsNew=list(set(EdgeLsNew))	
         print(Dec2AncNew, '\n',Anc2DecNew, '\n',Bra2MutCNew, '\n',EdgeLsNew)	
         return Dec2AncNew, Anc2DecNew, Bra2MutCNew, EdgeLsNew
     def GetNeighboringBranch(self,Bra,Bra2SNVnum, Dec2Anc, Anc2Dec):	
             Neigh=[] 		
             TotMut=Bra2SNVnum[Bra]
             NeighCombo=[]
             NeighComboU=[]		
             if (Bra in Dec2Anc)==True: #Anc and Sib Up
                 Anc=Dec2Anc[Bra]
                 TotMut+=Bra2SNVnum[Anc]	
     		
                 Neigh.append([Anc])
                 NeighCombo.append(Anc)			
                 Decs=Anc2Dec[Anc]
     		
                 for Dec in Decs:
                    if Dec!=Bra:
                        Neigh.append([Dec])
                        NeighCombo.append(Dec)				   
                        TotMut+=Bra2SNVnum[Dec]
                 if len(NeighCombo)>1: Neigh.append(NeighCombo)
             if (Bra in Anc2Dec)==True: #Dec Up  
                 Decs=Anc2Dec[Bra]
     
                 for Dec in Decs:
                        Neigh.append([Dec])
                        NeighComboU.append(Dec)				   
                        TotMut+=Bra2SNVnum[Dec]
                 if len(NeighComboU)>1: Neigh.append(NeighComboU)
             if len(NeighComboU)>1 and len(NeighCombo)>1: Neigh.append(NeighComboU+NeighCombo) 
             return Neigh, TotMut
     def MakeMutCountFile(self,BraCombo, Bra2File,OutFName):
         print(BraCombo)	
         print(Bra2File)		 
         BraLs=BraCombo.split('B')
         FileLs=[]
         for Bra in BraLs:
            if Bra!='':
               Bra='B'+Bra	   
               FileLs.append(Bra2File[Bra])
         print(FileLs,OutFName)			   
         self.pool(FileLs,OutFName)		  
     	
     def pool(self,FileLs,OutFileName):
            #import os
            #dir=os.getcwd()			
            SaveFile=file_writer()	 
            Sig2Count={}
            SigOrder=[]
            First='y'
            MutT=0    
            for i in FileLs:
                print(i)  
               # i=Fol+i+'.csv'
                i=open(i,'r').readlines()
                for Line in i:	   
                  Line=Line.strip().split(',')
                  if First=='y': 
                      SigOrder.append(Line[0])
                      Sig2Count[Line[0]]=0
                  Sig2Count[Line[0]]+=int(Line[1])
                  MutT+=int(Line[1])				  
                First='n'
            
            out=''
            for Sig in SigOrder:
                out+=Sig+','+str(Sig2Count[Sig])+'\n'
            print(MutT)				
            SaveFile.GetOut(OutFileName,out)	
     