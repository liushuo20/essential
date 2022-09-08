from Bio import SeqIO
import sys
if __name__ == '__main__':
    LenArgv=len(sys.argv)
    output='yn00.format' #default out name
    for i in range(LenArgv):
        if '-i'==sys.argv[i]:
            input=sys.argv[i+1]
        if '-o'==sys.argv[i]:
            output=sys.argv[i+1]
        if '-h'==sys.argv[i]:
            print 'help information:'
            print '-i   Your aligned file that you input'
            print '-o   Name of output file. Default yn00.format'
            exit()
    out=open(output, 'w')
    alinged=SeqIO.parse(open(input),'fasta')
    seqDic={}
    for j in alinged:
        seqDic[str(j.id)]=str(j.seq)[:-3]
    countNum=len(seqDic.keys())
    #check sequence length
    seqLenList=[]
    for k in seqDic.keys():
        seqLenList.append(len(seqDic[k]))
    if len(list(set(seqLenList))) != 1:
        print "The sequences in your input file were not aligned"
        exit()
    out.write(str(countNum)+"  "+str(seqLenList[0])+"\n")
    for m in seqDic.keys():
        out.write(m+"  "+seqDic[m]+"\n")
    out.close()
    print "Format transform was finished"








