"""Independent anchor/region numbering prototype, standard library only.

No ANARCI/ANARCII imports, neural weights, HTTP, or query-sequence lookup.
Public numbered frameworks supply alignment landmarks; loop labels are computed.
"""
from pathlib import Path
import json
import re
import string

SCHEMES = ('kabat', 'chothia', 'imgt', 'martin', 'aho')


class UnsupportedNumbering(ValueError):
    pass


def align(reference, query, deletion_penalties=None):
    """Needleman-Wunsch traceback written here; deterministic diagonal-first ties."""
    n, m = len(reference), len(query)
    deletion_penalties=deletion_penalties or [-5.0]*n
    score = [[0.0]*(m+1) for _ in range(n+1)]
    trace = [[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1):score[i][0]=score[i-1][0]+deletion_penalties[i-1];trace[i][0]=1
    for j in range(1,m+1):score[0][j]=-5*j;trace[0][j]=2
    for i in range(1,n+1):
        for j in range(1,m+1):
            a,b=reference[i-1],query[j-1]
            substitution=0 if 'X' in (a,b) else (3 if a==b else -1)
            options=(score[i-1][j-1]+substitution,score[i-1][j]+deletion_penalties[i-1],score[i][j-1]-5)
            best=max(range(3),key=lambda k:options[k]);score[i][j]=options[best];trace[i][j]=best
    pairs=[];i,j=n,m
    while i or j:
        k=trace[i][j]
        if k==0:pairs.append((i-1,j-1));i-=1;j-=1
        elif k==1:pairs.append((i-1,None));i-=1
        else:pairs.append((None,j-1));j-=1
    return score[n][m],list(reversed(pairs))


def centered(start,end,count):
    """IMGT symmetric gaps; excess labels use website decimal insertion notation."""
    width=end-start+1
    if count<=width:
        left=(count+1)//2;right=count//2
        return [str(x) for x in range(start,start+left)]+[str(x) for x in range(end-right+1,end+1)]
    if (start,end)!=(105,117):raise UnsupportedNumbering('Loop exceeds supported IMGT width')
    extra=count-width;left=extra//2;right=extra-left
    return ([str(x)for x in range(105,112)]+['111.'+str(x)for x in range(1,left+1)]
            +['112.'+str(x)for x in range(right,0,-1)]+[str(x)for x in range(112,118)])


def linear(start,end):return [str(x)for x in range(start,end+1)]


def loop_labels(chain,region,length,scheme,complete=True,light_subtype="kappa"):
    if chain=='H':
        if region==1:
            if not 10<=length<=12:raise UnsupportedNumbering('H1 supported lengths 10..12')
            extra=length-10
            if scheme=='imgt':return centered(27,38,length-2)+['39','40']
            if scheme=='aho':return linear(27,26+length-4)+linear(39,42)
            split=35 if scheme=='kabat' else 31
            return linear(26,split)+[str(split)+string.ascii_uppercase[i]for i in range(extra)]+linear(split+1,35)
        if region==2:
            if not 15<=length<=19:raise UnsupportedNumbering('H2 supported lengths 15..19')
            if scheme=='imgt':return ['55']+centered(56,65,length-9)+linear(66,72)+['74']
            if scheme=='aho':
                # Reproduced in a public PDB6VJA insertion control, not a query patch.
                # Legacy abYsis truncates this conversion at the 19-residue case.
                if length==19:return ['57','58','59','60','61','65']+['']*13
                if length==15:return linear(57,59)+linear(65,76)
                return linear(57,56+length-13)+linear(64,76)
            if length==15:return ['50','51','52']+linear(54,65)
            return ['50','51','52']+['52'+string.ascii_uppercase[i]for i in range(length-16)]+linear(53,65)
        if region==3:
            if scheme=='imgt':return centered(105,117,length)
            if scheme=='aho':
                if complete:
                    if length<2 or length>32:raise UnsupportedNumbering('Unsupported AHo H3 length')
                    return linear(107,106+length-2)+['137','138']
                return linear(107,106+length)
            if not complete:return linear(93,92+length)
            if length<5:raise UnsupportedNumbering('Complete H3 needs at least 3 loop residues')
            n=length-2
            if n<=8:return ['93','94']+linear(95,94+n-2)+['101','102']
            if n>34:raise UnsupportedNumbering('H3 insertion alphabet exhausted')
            return ['93','94']+linear(95,100)+['100'+string.ascii_uppercase[i]for i in range(n-8)]+['101','102']
    else:
        if region==1:
            if not 10<=length<=17:raise UnsupportedNumbering('L1 supported lengths 10..17')
            if scheme=='imgt':return ['24','25','26']+centered(27,38,length-5)+['39','40']
            if scheme=='aho':
                if light_subtype=='lambda':return linear(24,27)+linear(29,28+length-8)+linear(39,42)
                tail=5 if length>=15 else 4
                return linear(24,26)+linear(29,28+length-3-tail)+linear(43-tail,42)
            if length==10:
                deleted={'kabat':28,'chothia':31,'martin':30}[scheme]
                return [str(x)for x in range(24,35)if x!=deleted]
            site=27 if scheme=='kabat' else 30
            return linear(24,site)+[str(site)+string.ascii_uppercase[i]for i in range(length-11)]+linear(site+1,34)
        if region==2:
            if length!=7:raise UnsupportedNumbering('Noncanonical L2 not yet supported')
            if scheme=='imgt':return ['56','65','66','67','68','69','70']
            if scheme=='aho':return ['58','59','60','69','70','71','72']
            return linear(50,56)
        if region==3:
            if scheme=='imgt':return centered(105,117,length)
            if scheme=='aho':return linear(107,106+length-2)+['137','138'] if complete else linear(107,106+length)
            if not complete:return linear(89,88+length)
            if not 7<=length<=15:raise UnsupportedNumbering('L3 supported lengths 7..15')
            if length<=9:return linear(89,88+length-2)+['96','97']
            return linear(89,95)+['95'+string.ascii_uppercase[i]for i in range(length-9)]+['96','97']
    raise UnsupportedNumbering('Unknown loop')


class AnchorNumberer:
    def __init__(self, framework_asset):
        data=json.loads(Path(framework_asset).read_text())
        if data['schema']!='public_numbered_frameworks_v2':raise ValueError('Unsupported public framework schema')
        self.templates=data['templates']
        assert {t['chain'] for t in self.templates.values()}=={'H','L'}
        for _,t in self.templates.items():
            chain=t['chain']
            assert set(t['labels'])==set(SCHEMES)
            assert all(len(ps)==len(t['sequence']) for ps in t['labels'].values())
            for x,y in t['segments'].values():
                assert 0<=x<y<=len(t['sequence'])
                assert all(t['labels'][sc][i].startswith(chain) for sc in SCHEMES for i in range(x,y))

    def framework_alignment(self,t,key,query):
        x,y=t['segments'][key];ref=t['sequence'][x:y]
        deletion=[-5.0]*len(ref)
        if key=='fr1' and len(query)<len(ref):
            # N-terminal truncations are not internal framework deletions.
            for i in range(len(ref)-len(query)):deletion[i]=-.1
        if key=='fr2' and t['chain']=='H':
            signature=tuple(t['labels']['kabat'][x:y])
            refs=[]
            for other in self.templates.values():
                u,v=other['segments'][key]
                if other['chain']=='H' and tuple(other['labels']['kabat'][u:v])==signature:
                    refs.append(align(other['sequence'][u:v],query,deletion))
            return max(refs,key=lambda item:item[0])
        if key=='fr4' and len(query)<len(ref):
            # A C-terminal truncation cannot right-align itself to terminal labels.
            return 0.0,[(i,i) for i in range(len(query))]
        return align(ref,query,deletion)

    def locate(self,sequence):
        choices=[]
        terminals=[i for i in range(2,len(sequence)) if sum(a==b for a,b in zip(sequence[i-2:i+1],'YYC'))>=2]
        for template_id,t in self.templates.items():
            chain=t['chain'];cached={}
            def framework_score(key,a,b):
                cachekey=(key,a,b)
                if cachekey not in cached:cached[cachekey]=self.framework_alignment(t,key,sequence[a:b])[0]
                return cached[cachekey]
            for c1 in range(min(71,len(sequence))):
                start=max(0,c1-t['c1']);end1=c1+(4 if chain=='H' else 1)
                for loop1length in (range(10,13) if chain=='H' else range(10,18)):
                    w1=end1+loop1length
                    for fr2length in ([13,14] if chain=='H' else [15]):
                        afterfr2=w1+fr2length
                        for c2 in terminals:
                            for fr3length in ([30] if chain=='H' else [32]):
                                startfr3=c2+1-fr3length;n2=startfr3-afterfr2
                                if not (15<=n2<=19 if chain=='H' else n2==7):continue
                                if c2+1>=len(sequence):continue
                                segments={'fr1':(start,end1),'fr2':(w1,afterfr2),'fr3':(startfr3,c2+1)}
                                if any(a<0 or b>len(sequence) or b<=a for a,b in segments.values()):continue
                                scores=[framework_score(key,a,b) for key,(a,b) in segments.items()]
                                lengths=[b-a for a,b in segments.values()]
                                if any(sc/n<.2 for sc,n in zip(scores,lengths)):continue
                                anchor_support=(6 if sequence[c1]=='C' else 0)+(3 if sequence[w1]=='W' else 0)
                                subtype_support=0
                                if chain=='L':
                                    tail=re.search(r'[FW]G[A-Z]G',sequence[c2+1:])
                                    if tail:
                                        start4=c2+1+tail.start();fr4seq=sequence[start4:start4+10]
                                        lambda_pattern=bool(re.fullmatch(r'[FW]G.GT.[LIV][TS]V[LV]',fr4seq))
                                        kappa_pattern=bool(re.fullmatch(r'[FW]G.GT..[ED][ILV][KR]',fr4seq))
                                        is_lambda=t.get('light_subtype')=='lambda'
                                        if lambda_pattern:subtype_support=12 if is_lambda else -12
                                        elif kappa_pattern:subtype_support=-12 if is_lambda else 12
                                score=(sum(scores)+anchor_support+subtype_support)/sum(lengths)
                                if sum(scores)/sum(lengths)<.45:continue
                                # Prefer the complete framework when support is tied.
                                choices.append((score,chain,segments,c1,c2,w1,template_id))
        if not choices:return None
        return max(choices,key=lambda x:(x[0],sum(b-a for a,b in x[2].values()),x[1],x[6]))

    def number(self,sequence):
        sequence=sequence.upper().replace('-','').replace('?','X')
        if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWYX'):raise ValueError('Invalid amino acid input')
        found=self.locate(sequence)
        if found is None:return {'status':'unclassified','schemes':{s:['']*len(sequence)for s in SCHEMES}}
        score,chain,segments,c1,c2,w1,template_id=found;t=self.templates[template_id]
        h3start=c2+1
        motif=r'W[GAX][A-Z]G' if chain=='H' else r'[FW]G[A-Z]G'
        match=re.search(motif,sequence[h3start:]);fr4=h3start+match.start() if match else None
        complete=fr4 is not None
        loop3end=fr4 if complete else len(sequence)
        if complete:
            width=t['segments']['fr4'][1]-t['segments']['fr4'][0]
            segments['fr4']=(fr4,min(len(sequence),fr4+width))
        loops={1:(segments['fr1'][1],w1),2:(segments['fr2'][1],segments['fr3'][0]),3:(h3start,loop3end)}
        result={};notes=[]
        for scheme in SCHEMES:
            positions=['']*len(sequence)
            try:
                for key,(a,b) in segments.items():
                    x,y=t['segments'][key];_,pairs=self.framework_alignment(t,key,sequence[a:b])
                    for i,j in pairs:
                        if i is None and j is not None:raise UnsupportedNumbering('Framework insertion outside validated anchors')
                        if i is not None and j is not None:positions[a+j]=t['labels'][scheme][x+i]
                for region,(a,b) in loops.items():
                    if chain=='L' and region==1 and scheme=='imgt':
                        # IMGT L24..26 precedes the centered L27..38 loop; L39/40 follows.
                        labels=['24','25','26']+centered(27,38,b-a-5)+['39','40']
                    else:labels=loop_labels(chain,region,b-a,scheme,complete=complete if region==3 else True,light_subtype=t.get('light_subtype','kappa'))
                    assert len(labels)==b-a,(chain,scheme,region,len(labels),b-a)
                    for i,label in enumerate(labels,a):positions[i]=chain+label if label else ''
                if chain=='H' and scheme=='aho' and loops[2][1]-loops[2][0]==19:
                    positions[loops[2][0]+6:]=['']*(len(sequence)-loops[2][0]-6)
                    notes.append('Legacy abYsis AHo partial output for H2 length19; confirmed on public insertion control')
            except UnsupportedNumbering as error:
                result[scheme]={'status':'unsupported','reason':str(error),'positions':['']*len(sequence)};continue
            result[scheme]={'status':'numbered','positions':positions}
        return {'status':'numbered','chain':chain,'framework_template':template_id,'alignment_score_not_probability':score,'segments':segments,'loops':loops,
                'complete_fr4':complete,'schemes':result,'compatibility_notes':notes}
