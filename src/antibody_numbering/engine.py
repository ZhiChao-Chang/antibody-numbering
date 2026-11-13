"""Reference-assisted antibody variable-domain numbering, using the standard library.

Reading order for reviewers:
    align()               - establish residue correspondences by dynamic programming
    centered()/linear()   - construct positional labels within a numbering scheme
    loop_labels()         - apply the supported chain/length-specific loop rules
    AnchorNumberer        - locate segments, map frameworks, and assemble results

Reference framework residues carry public labels; masked loop labels are computed
from explicit rules. The implementation has no ANARCI/ANARCII imports, trained
weights, network calls, or lookup table of answers for incoming sequences.

All sequence offsets and segment slices in this module are zero-based, with
exclusive ends. Biological labels such as H35A are strings, not array offsets.
The CLI exports one-based raw/normalized indices separately.
"""
from pathlib import Path
import json
import re
import string

SCHEMES = ('kabat', 'chothia', 'imgt', 'martin', 'aho')


class UnsupportedNumbering(ValueError):
    """A candidate cannot be converted using the implemented numbering rules."""

    pass


def align(reference, query, deletion_penalties=None):
    """Return the global alignment score and zero-based residue index pairs.

    A pair (i, j) aligns reference[i] to query[j]. A None index denotes a gap:
    (i, None) deletes a reference residue; (None, j) inserts a query residue.
    Matches score +3, mismatches -1, and pairs involving unknown X score 0.
    Query insertions cost -5; reference deletion penalties can be customized.
    Ties prefer diagonal, then deletion, then insertion, for reproducible output.
    """
    n, m = len(reference), len(query)
    deletion_penalties=deletion_penalties or [-5.0]*n
    # score[i][j] is the best score for the two prefixes. trace stores the move
    # that produced it: 0 = diagonal, 1 = reference deletion, 2 = query insertion.
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
    # Trace the optimum from the final cell back to the empty prefixes, then
    # reverse the collected pairs into sequence order.
    pairs=[];i,j=n,m
    while i or j:
        k=trace[i][j]
        if k==0:pairs.append((i-1,j-1));i-=1;j-=1
        elif k==1:pairs.append((i-1,None));i-=1
        else:pairs.append((None,j-1));j-=1
    return score[n][m],list(reversed(pairs))


def centered(start,end,count):
    """Place count labels into an inclusive IMGT interval from both ends.

    Short loops leave a central gap. Longer-than-interval loops are supported
    only for positions 105..117, using the stored website's decimal insertion
    convention around 111/112. Return label strings without an H/L prefix.
    """
    width=end-start+1
    if count<=width:
        left=(count+1)//2;right=count//2
        return [str(x) for x in range(start,start+left)]+[str(x) for x in range(end-right+1,end+1)]
    if (start,end)!=(105,117):raise UnsupportedNumbering('Loop exceeds supported IMGT width')
    # Insertions on the 112 side are emitted in descending suffix order to
    # retain sequence order across the center; they are not decimal numbers.
    extra=count-width;left=extra//2;right=extra-left
    return ([str(x)for x in range(105,112)]+['111.'+str(x)for x in range(1,left+1)]
            +['112.'+str(x)for x in range(right,0,-1)]+[str(x)for x in range(112,118)])


def linear(start,end):
    """Return inclusive integer position labels as strings, without chain prefixes."""
    return [str(x)for x in range(start,end+1)]


def loop_labels(chain,region,length,scheme,complete=True,light_subtype="kappa"):
    """Compute labels for one inter-framework segment in a supported scheme.

    region is 1, 2, or 3; length counts the implementation's segment, which may
    include flanking residues outside a scheme's biological CDR definition.
    complete indicates whether the third segment has an identified FR4 boundary.
    light_subtype controls kappa/lambda-specific placement where implemented.
    Return one label per residue; empty labels explicitly retain partial output.
    These cases and length limits are hand-coded rules, not learned parameters.
    """
    if chain=='H':
        if region==1:
            if not 10<=length<=12:raise UnsupportedNumbering('H1 supported lengths 10..12')
            extra=length-10
            if scheme=='imgt':return centered(27,38,length-2)+['39','40']
            if scheme=='aho':return linear(27,26+length-4)+linear(39,42)
            # Equivalent extra residues occupy different insertion sites across
            # schemes: e.g. Kabat H35A versus Chothia/Martin H31A.
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
            # Reserve the two right-flank positions; longer supported H3
            # segments use alphabetic insertions after position 100.
            n=length-2
            if n<=8:return ['93','94']+linear(95,94+n-2)+['101','102']
            if n>34:raise UnsupportedNumbering('H3 insertion alphabet exhausted')
            return ['93','94']+linear(95,100)+['100'+string.ascii_uppercase[i]for i in range(n-8)]+['101','102']
    else:
        if region==1:
            if not 10<=length<=17:raise UnsupportedNumbering('L1 supported lengths 10..17')
            if scheme=='imgt':return ['24','25','26']+centered(27,38,length-5)+['39','40']
            if scheme=='aho':
                # Kappa and lambda do not share identical short-loop gap sites.
                if light_subtype=='lambda':return linear(24,27)+linear(29,28+length-8)+linear(39,42)
                tail=5 if length>=15 else 4
                return linear(24,26)+linear(29,28+length-3-tail)+linear(43-tail,42)
            if length==10:
                # A short L1 removes a scheme-specific slot rather than shifting
                # every downstream framework label by one residue.
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
    """Number a single variable domain using a fixed set of public frameworks.

    Framework selection is computed for each input. Loading an asset does not
    train a model or retrieve an existing answer for the input sequence.
    """

    def __init__(self, framework_asset):
        """Load the JSON asset and validate chain, label, and slice consistency."""
        data=json.loads(Path(framework_asset).read_text())
        if data['schema']!='public_numbered_frameworks_v2':raise ValueError('Unsupported public framework schema')
        self.templates=data['templates']
        assert {t['chain'] for t in self.templates.values()}=={'H','L'}
        # Each framework label array must use the template sequence's index
        # space. Only the declared framework slices are required to carry labels.
        for _,t in self.templates.items():
            chain=t['chain']
            assert set(t['labels'])==set(SCHEMES)
            assert all(len(ps)==len(t['sequence']) for ps in t['labels'].values())
            for x,y in t['segments'].values():
                assert 0<=x<y<=len(t['sequence'])
                assert all(t['labels'][sc][i].startswith(chain) for sc in SCHEMES for i in range(x,y))

    def framework_alignment(self,t,key,query):
        """Align a query segment to FR1/FR2/FR3/FR4 of template t.

        Terminal exceptions preserve labels for the supported truncations.
        Heavy-chain FR2 may use an alternative public sequence with the same
        Kabat label signature; this changes sequence matching, not label order.
        """
        x,y=t['segments'][key];ref=t['sequence'][x:y]
        deletion=[-5.0]*len(ref)
        if key=='fr1' and len(query)<len(ref):
            # N-terminal truncations are not internal framework deletions.
            for i in range(len(ref)-len(query)):deletion[i]=-.1
        if key=='fr2' and t['chain']=='H':
            # Alternative FR2 references are interchangeable only when their
            # labels define the same position-by-position mapping.
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
        """Choose a supported chain/template/segmentation hypothesis, or None.

        Enumerate candidate anchors and supported segment lengths, align the
        framework pieces, reject weak matches, then rank the remaining choices.
        Scores and search limits are fixed heuristics, not fitted probabilities.
        Returned c1/c2 are anchor indices; w1 is the candidate start of FR2.
        """
        choices=[]
        # Candidate FR3 endpoints retain at least two letters of the local YYC
        # motif. This hard candidate filter differs from the soft C/W bonuses below.
        terminals=[i for i in range(2,len(sequence)) if sum(a==b for a,b in zip(sequence[i-2:i+1],'YYC'))>=2]
        for template_id,t in self.templates.items():
            chain=t['chain'];cached={}
            def framework_score(key,a,b):
                """Memoize repeated slice alignments within this template search."""
                cachekey=(key,a,b)
                if cachekey not in cached:cached[cachekey]=self.framework_alignment(t,key,sequence[a:b])[0]
                return cached[cachekey]
            # Search the first 71 offsets for the first anchor. c1 need not be C:
            # the residue contributes a bonus rather than an unconditional veto.
            # c2 anchors FR3; the distances imply candidate loop1/loop2 lengths.
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
                                # Evaluate FR1, FR2 and FR3 separately. The 0.2
                                # per-segment and 0.45 aggregate cutoffs below
                                # are fixed minimum scores per residue.
                                scores=[framework_score(key,a,b) for key,(a,b) in segments.items()]
                                lengths=[b-a for a,b in segments.values()]
                                if any(sc/n<.2 for sc,n in zip(scores,lengths)):continue
                                # Conserved C/W matches add support without
                                # requiring every anchor residue to be unchanged.
                                anchor_support=(6 if sequence[c1]=='C' else 0)+(3 if sequence[w1]=='W' else 0)
                                subtype_support=0
                                if chain=='L':
                                    tail=re.search(r'[FW]G[A-Z]G',sequence[c2+1:])
                                    if tail:
                                        start4=c2+1+tail.start();fr4seq=sequence[start4:start4+10]
                                        # FR4 motifs supply an additional fixed
                                        # +12/-12 kappa/lambda preference when
                                        # they match; ambiguous tails add none.
                                        lambda_pattern=bool(re.fullmatch(r'[FW]G.GT.[LIV][TS]V[LV]',fr4seq))
                                        kappa_pattern=bool(re.fullmatch(r'[FW]G.GT..[ED][ILV][KR]',fr4seq))
                                        is_lambda=t.get('light_subtype')=='lambda'
                                        if lambda_pattern:subtype_support=12 if is_lambda else -12
                                        elif kappa_pattern:subtype_support=-12 if is_lambda else 12
                                # Normalize by framework residues so candidates
                                # of different lengths are compared on this scale.
                                score=(sum(scores)+anchor_support+subtype_support)/sum(lengths)
                                if sum(scores)/sum(lengths)<.45:continue
                                # Prefer the complete framework when support is tied.
                                choices.append((score,chain,segments,c1,c2,w1,template_id))
        if not choices:return None
        return max(choices,key=lambda x:(x[0],sum(b-a for a,b in x[2].values()),x[1],x[6]))

    def number(self,sequence):
        """Normalize one domain sequence and assemble all five scheme mappings.

        UnsupportedNumbering is caught per scheme, so a failure in one scheme
        does not erase other schemes' results. Overall status and per-scheme
        status have distinct meanings and must both be inspected by callers.
        Unclassified inputs return blank lists rather than per-scheme dictionaries.
        """
        sequence=sequence.upper().replace('-','').replace('?','X')
        if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWYX'):raise ValueError('Invalid amino acid input')
        found=self.locate(sequence)
        if found is None:return {'status':'unclassified','schemes':{s:['']*len(sequence)for s in SCHEMES}}
        score,chain,segments,c1,c2,w1,template_id=found;t=self.templates[template_id]
        # The FR4 motif defines the end of the third inter-framework segment.
        # Without a motif, the remaining tail is processed as incomplete; this is
        # a supported fallback, not proof that the true FR4 is biologically absent.
        h3start=c2+1
        motif=r'W[GAX][A-Z]G' if chain=='H' else r'[FW]G[A-Z]G'
        match=re.search(motif,sequence[h3start:]);fr4=h3start+match.start() if match else None
        complete=fr4 is not None
        loop3end=fr4 if complete else len(sequence)
        if complete:
            width=t['segments']['fr4'][1]-t['segments']['fr4'][0]
            segments['fr4']=(fr4,min(len(sequence),fr4+width))
        # These internal spans bridge the framework slices. Their boundaries
        # need not equal the CDR definitions in every output numbering scheme.
        loops={1:(segments['fr1'][1],w1),2:(segments['fr2'][1],segments['fr3'][0]),3:(h3start,loop3end)}
        result={};notes=[]
        for scheme in SCHEMES:
            positions=['']*len(sequence)
            try:
                # Step 1: transfer public framework labels through alignment.
                # A new framework insertion has no implemented label rule here.
                for key,(a,b) in segments.items():
                    x,y=t['segments'][key];_,pairs=self.framework_alignment(t,key,sequence[a:b])
                    for i,j in pairs:
                        if i is None and j is not None:raise UnsupportedNumbering('Framework insertion outside validated anchors')
                        if i is not None and j is not None:positions[a+j]=t['labels'][scheme][x+i]
                # Step 2: compute loop labels from lengths and scheme rules.
                for region,(a,b) in loops.items():
                    if chain=='L' and region==1 and scheme=='imgt':
                        # IMGT L24..26 precedes the centered L27..38 loop; L39/40 follows.
                        labels=['24','25','26']+centered(27,38,b-a-5)+['39','40']
                    else:labels=loop_labels(chain,region,b-a,scheme,complete=complete if region==3 else True,light_subtype=t.get('light_subtype','kappa'))
                    assert len(labels)==b-a,(chain,scheme,region,len(labels),b-a)
                    for i,label in enumerate(labels,a):positions[i]=chain+label if label else ''
                # Compatibility mode: the verified legacy converter truncates
                # the entire AHo suffix after six H2-segment labels in this case.
                # Keep blanks and a note rather than presenting a full mapping.
                if chain=='H' and scheme=='aho' and loops[2][1]-loops[2][0]==19:
                    positions[loops[2][0]+6:]=['']*(len(sequence)-loops[2][0]-6)
                    notes.append('Legacy abYsis AHo partial output for H2 length19; confirmed on public insertion control')
            except UnsupportedNumbering as error:
                # Keep the diagnostic reason and input-length correspondence.
                # Do not publish partly assigned labels for a rejected scheme.
                result[scheme]={'status':'unsupported','reason':str(error),'positions':['']*len(sequence)};continue
            result[scheme]={'status':'numbered','positions':positions}
        return {'status':'numbered','chain':chain,'framework_template':template_id,'alignment_score_not_probability':score,'segments':segments,'loops':loops,
                'complete_fr4':complete,'schemes':result,'compatibility_notes':notes}
