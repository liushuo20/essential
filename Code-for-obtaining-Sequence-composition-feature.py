# -*- coding: cp936 -*-

########################################################################
#def usage():
def subSeq(seq,w,k):
    firstPhaseSubSeq=[]
    secondPhaseSubSeq=[]
    thirdPhaseSubSeq=[]
    PhaseSubSeq=[]

    if w!=1:
        firstPhaseIndex=range(0,len(seq),3)
        #print firstPhaseIndex
        secondPhaseIndex=range(1,len(seq),3)
        thirdPhaseIndex=range(2,len(seq),3)
        seqIndex=range(0,len(seq))
        for i1 in firstPhaseIndex:
            if w-2+i1+k in seqIndex:
                #print w-2+i1+k
                firstPhaseSubSeq.append(seq[i1:(i1+w-1)]+seq[w-2+i1+k])
        for i2 in secondPhaseIndex:
            if w-2+i2+k in seqIndex:
                
                secondPhaseSubSeq.append(seq[i2:(i2+w-1)]+seq[w-2+i2+k])
        for i3 in thirdPhaseIndex:
            #print thirdPhaseIndex
            if w-2+i3+k in seqIndex:
                thirdPhaseSubSeq.append(seq[i3:(i3+w-1)]+seq[w-2+i3+k])
                #print 'sss'
        PhaseSubSeq.append(firstPhaseSubSeq)
        PhaseSubSeq.append(secondPhaseSubSeq)
        PhaseSubSeq.append(thirdPhaseSubSeq)
        #print PhaseSubSeq
        return PhaseSubSeq
        
    else:
        firstPhaseIndex=range(0,len(seq),3)
        secondPhaseIndex=range(1,len(seq),3)
        thirdPhaseIndex=range(2,len(seq),3)
        for i1 in firstPhaseIndex:
            firstPhaseSubSeq.append(seq[i1])
        for i2 in secondPhaseIndex:
            secondPhaseSubSeq.append(seq[i2])
        for i3 in thirdPhaseIndex:
            thirdPhaseSubSeq.append(seq[i3])
        PhaseSubSeq.append(firstPhaseSubSeq)
        PhaseSubSeq.append(secondPhaseSubSeq)
        PhaseSubSeq.append(thirdPhaseSubSeq)
        return PhaseSubSeq
    
    

def Nucleiotide2x(baseList,w):
    #print w
    base=['A','T','G','C']
    xValue=[]
    if w==1:
        num1=float(baseList.count('A'))+float(baseList.count('G'))
        num2=float(baseList.count('C'))+float(baseList.count('T'))
        xValue.append((num1-num2)/len(baseList))
        return xValue
    if w==2:
        #print 'ssss'
        for i in base:
            num1=float(baseList.count(i+'A'))+float(baseList.count(i+'G'))
            num2=float(baseList.count(i+'C'))+float(baseList.count(i+'T'))
            #print num2
            #print i+'A', i+'G', num1
            #print i+'C', i+'T', num2
            xValue.append((num1-num2)/len(baseList))
            #print xValue
        return xValue
    if w==3:
        for i in base:
            for j in base:
                num1=float(baseList.count(i+j+'A'))+float(baseList.count(i+j+'G'))
                num2=float(baseList.count(i+j+'C'))+float(baseList.count(i+j+'T'))
                xValue.append((num1-num2)/len(baseList))
        return xValue
    if w==4:
        for i in base:
            for j in base:
                for m in base:
                    num1=float(baseList.count(i+j+m+'A'))+float(baseList.count(i+j+m+'G'))
                    num2=float(baseList.count(i+j+m+'C'))+float(baseList.count(i+j+m+'T'))
                    xValue.append((num1-num2)/len(baseList))
        return xValue

def Nucleiotide2y(baseList,w):
    base=['A','T','G','C']
    yValue=[]
    if w==1:
        num1=float(baseList.count('A'))+float(baseList.count('C'))
        num2=float(baseList.count('G'))+float(baseList.count('T'))
        yValue.append((num1-num2)/len(baseList))
        return yValue
    if w==2:
        for i in base:
            num1=float(baseList.count(i+'A'))+float(baseList.count(i+'C'))
            num2=float(baseList.count(i+'G'))+float(baseList.count(i+'T'))
            yValue.append((num1-num2)/len(baseList))
            #print
        return yValue
    if w==3:
        for i in base:
            for j in base:
                num1=float(baseList.count(i+j+'A'))+float(baseList.count(i+j+'C'))
                num2=float(baseList.count(i+j+'G'))+float(baseList.count(i+j+'T'))
                yValue.append((num1-num2)/len(baseList))
        return yValue
    if w==4:
        for i in base:
            for j in base:
                for m in base:
                    num1=float(baseList.count(i+j+m+'A'))+float(baseList.count(i+j+m+'C'))
                    num2=float(baseList.count(i+j+m+'G'))+float(baseList.count(i+j+m+'T'))
                    yValue.append((num1-num2)/len(baseList))
        return yValue 
            

def Nucleiotide2z(baseList,w):
    base=['A','T','G','C']
    zValue=[]
    if w==1:
        num1=float(baseList.count('A'))+float(baseList.count('T'))
        num2=float(baseList.count('G'))+float(baseList.count('C'))
        zValue.append((num1-num2)/len(baseList))
        return zValue
    if w==2:
        for i in base:
            num1=float(baseList.count(i+'A'))+float(baseList.count(i+'T'))
            num2=float(baseList.count(i+'G'))+float(baseList.count(i+'C'))
            zValue.append((num1-num2)/len(baseList))
        return zValue
    if w==3:
        for i in base:
            for j in base:
                num1=float(baseList.count(i+j+'A'))+float(baseList.count(i+j+'T'))
                num2=float(baseList.count(i+j+'G'))+float(baseList.count(i+j+'C'))
                zValue.append((num1-num2)/len(baseList))
        return zValue
    if w==4:
        for i in base:
            for j in base:
                for m in base:
                    num1=float(baseList.count(i+j+m+'A'))+float(baseList.count(i+j+m+'T'))
                    num2=float(baseList.count(i+j+m+'G'))+float(baseList.count(i+j+m+'C'))
                    zValue.append((num1-num2)/len(baseList))
        return zValue         
        
  
def kwZcurve(seq,w,k):
    PhaseSubSeq=subSeq(seq,w,k)
    xf=[]
    yf=[]
    zf=[]
    seq2feature=[]
    for subPhase in PhaseSubSeq:
        x=Nucleiotide2x(subPhase,w)
        y=Nucleiotide2y(subPhase,w)
        z=Nucleiotide2z(subPhase,w)
        xf.extend(x)
        yf.extend(y)
        zf.extend(z)
    seq2feature.extend(xf)
    seq2feature.extend(yf)
    seq2feature.extend(zf)
    return seq2feature
    #print seq2feature

##############
#main program#
##############
if __name__=='__main__':
    import sys
    from Bio import SeqIO
    '''LenArgv=len(sys.argv)
    for i in range(LenArgv):
        if sys.argv[i]=='-w':
            windowSize=sys.argv[i+1]
        if sys.argv[i]=='-k':
            kGap=sys.argv[i+1]
        if sys.argv[i]=='-s':
            seq=sys.argv[i+1]'''
    
for endw in range (4,5):
    for endk in range(6,7):
        Variable=open('mus'+"w"+str(endw)+"_k"+str(endk),"w")
        Target_Record=SeqIO.parse(open('.\\testessential.txt','r'),'fasta')
        label='1'
        for record in Target_Record:
            feature=[]
            seqFurther=kwZcurve(str(record.seq),1,1)
            seqFurther=map(str,seqFurther)
            feature.extend(seqFurther)
            Variable.write(label+" ")
            for w in range(2,endw+1):
                for k in range(1,endk+1):
                
                    seqFurther=kwZcurve(str(record.seq),w,k)
                    seqFurther=map(str,seqFurther)
                    feature.extend(seqFurther)
            featureIndex=1
            for featureVar in feature:
                Variable.write(str(featureIndex)+":"+featureVar+" ")
                featureIndex=featureIndex+1
            Variable.write('\n')
    #===========================================================================================

        Target_Record=SeqIO.parse(open('.\\testessential.txt','r'),'fasta')
        label='-1'
        for record in Target_Record:
            feature=[]
            seqFurther=kwZcurve(str(record.seq),1,1)
            seqFurther=map(str,seqFurther)
            feature.extend(seqFurther)
            Variable.write(label+" ")
            for w in range(2,endw+1):
                for k in range(1,endk+1):
                    #print "dd"
                
                    seqFurther=kwZcurve(str(record.seq),w,k)
                    seqFurther=map(str,seqFurther)
                    feature.extend(seqFurther)
            featureIndex=1
            for featureVar in feature:
                Variable.write(str(featureIndex)+":"+featureVar+" ")
                featureIndex=featureIndex+1
            Variable.write('\n')
        Variable.close()
    
    
    
        
    
        
    
    
        
    
            
    
            
        
    
        
        
        
        
    
    
    
