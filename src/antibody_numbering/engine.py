"""Reference-assisted antibody variable-domain numbering, using the standard library.

Reading order for reviewers:
    align()               - establish residue correspondences by dynamic programming
    centered()/linear()   - construct positional labels within a numbering scheme
    loop_labels()         - apply the supported chain/length-specific loop rules
    AnchorNumberer        - locate segments, map frameworks, and assemble results

Reference framework residues carry public labels; loop labels are computed
from explicit chain-, length-, and scheme-specific rules.

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


AA_CODES={aa:i for i,aa in enumerate('XACDEFGHIKLMNPQRSTVWY')}


def pack_residues(sequence):
    """Encode five-bit residue fields and a separate unknown-residue mask.

    This is a compact representation for scoring short framework segments.
    A one-bit marker selects the least-significant bit of every residue field.
    """
    word=0;unknown=0
    for i,aa in enumerate(sequence):
        shift=5*i;word|=AA_CODES[aa]<<shift
        if aa=='X':unknown|=1<<shift
    markers=((1<<(5*len(sequence)))-1)//31
    return word,unknown,markers


def packed_substitution_score(reference,query):
    """Exactly reproduce +3 match / -1 mismatch / 0 unknown for equal lengths.

    Folding the five XOR bits into each field's lowest bit marks mismatches.
    Python's native integer bit_count then counts known matches/mismatches in C.
    The result is identical to explicit residue-by-residue scoring.
    """
    ref,ref_unknown,markers=reference;value,query_unknown,other_markers=query
    if markers!=other_markers:raise ValueError('Packed scores require equal lengths')
    delta=ref^value
    mismatches=(delta|(delta>>1)|(delta>>2)|(delta>>3)|(delta>>4))&markers
    known=markers&~(ref_unknown|query_unknown)
    return 3*(known&~mismatches).bit_count()-(known&mismatches).bit_count()


def ordered_framework_alignment(reference, query, deletion_penalties=None):
    """Align a framework while preserving a label for every query residue.

    Supported framework profiles may omit reference positions, but cannot map
    an inserted query letter without a numbered slot. Equal lengths therefore
    have a direct correspondence. Shorter queries use a narrow deletion-only DP;
    truly longer queries retain the general aligner's insertion evidence.
    """
    n,m=len(reference),len(query)
    deletion_penalties=deletion_penalties or [-5.0]*n
    if m>n:return align(reference,query,deletion_penalties)
    def substitution(a,b):return 0 if 'X'in(a,b)else 3 if a==b else -1
    if n==m:return sum(substitution(a,b)for a,b in zip(reference,query)),[(i,i)for i in range(m)]
    gap=n-m;negative=-1e12
    previous=[negative]*(gap+1);previous[0]=0.0
    trace=[[0]*(gap+1)for _ in range(n+1)]
    for i in range(1,n+1):
        current=[negative]*(gap+1)
        for deleted in range(max(0,i-m),min(gap,i)+1):
            j=i-deleted
            diagonal=previous[deleted]+substitution(reference[i-1],query[j-1])if j else negative
            deletion=previous[deleted-1]+deletion_penalties[i-1]if deleted else negative
            if diagonal>=deletion:current[deleted]=diagonal;trace[i][deleted]=0
            else:current[deleted]=deletion;trace[i][deleted]=1
        previous=current
    pairs=[];i,j=n,m;deleted=gap
    while i:
        if trace[i][deleted]==0:pairs.append((i-1,j-1));i-=1;j-=1
        else:pairs.append((i-1,None));i-=1;deleted-=1
    return previous[gap],list(reversed(pairs))


def centered(start,end,count):
    """Place count labels into an inclusive IMGT interval from both ends.

    Short loops leave a central gap. Longer-than-interval loops are supported
    for H1/H2 intervals and for positions 105..117. The third interval uses
    the stored website's insertion convention around 111/112. Return strings
    without an H/L prefix.
    """
    width=end-start+1
    if count<=width:
        left=(count+1)//2;right=count//2
        return [str(x) for x in range(start,start+left)]+[str(x) for x in range(end-right+1,end+1)]
    if (start,end) in [(27,38),(56,65)]:
        midpoint=(start+end-1)//2
        return linear(start,midpoint)+[str(midpoint)+'.'+str(i)for i in range(1,count-width+1)]+linear(midpoint+1,end)
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
    Supported cases and length limits are defined explicitly below.
    """
    if chain=='H':
        if region==1:
            if not 7<=length<=16:raise UnsupportedNumbering('H1 supported lengths 7..16')
            extra=length-10
            if length<10 and scheme in ('kabat','chothia','martin'):
                return linear(26,25+length)if scheme=='kabat'else linear(26,21+length)+linear(32,35)
            if scheme=='imgt':return centered(27,38,length-2)+['39','40']
            if scheme=='aho':return linear(27,26+length-4)+linear(39,42)
            # Equivalent extra residues occupy different insertion sites across
            # schemes: e.g. Kabat H35A versus Chothia/Martin H31A.
            split=35 if scheme=='kabat' else 31
            return linear(26,split)+[str(split)+string.ascii_uppercase[i]for i in range(extra)]+linear(split+1,35)
        if region==2:
            if not 15<=length<=24:raise UnsupportedNumbering('H2 supported lengths 15..24')
            if scheme=='imgt':return ['55']+centered(56,65,length-9)+linear(66,72)+['74']
            if scheme=='aho':
                # Retain the supported legacy AHo partial-conversion rule.
                # Verified legacy abYsis truncation applies to supported lengths >=19.
                if length>=19:return ['57','58','59','60','61','65']+['']*(length-6)
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
            if not complete:
                if length<=8:return linear(93,92+length)
                if length>34:raise UnsupportedNumbering('Incomplete H3 insertion alphabet exhausted')
                return linear(93,100)+['100'+string.ascii_uppercase[i]for i in range(length-8)]
            if length<5:raise UnsupportedNumbering('Complete H3 needs at least 3 loop residues')
            # Reserve the two right-flank positions; longer supported H3
            # segments use alphabetic insertions after position 100.
            n=length-2
            if n<=8:return ['93','94']+linear(95,94+n-2)+['101','102']
            if n>34:raise UnsupportedNumbering('H3 insertion alphabet exhausted')
            return ['93','94']+linear(95,100)+['100'+string.ascii_uppercase[i]for i in range(n-8)]+['101','102']
    else:
        if region==1:
            if not 8<=length<=20:raise UnsupportedNumbering('L1 supported lengths 8..20')
            if scheme=='imgt':return ['24','25','26']+centered(27,38,length-5)+['39','40']
            if scheme=='aho':
                if length>(18 if light_subtype=='lambda'else 17):raise UnsupportedNumbering('AHo L1 exceeds validated slots')
                # Kappa and lambda do not share identical short-loop gap sites.
                if light_subtype=='lambda'and length==11:return linear(24,26)+linear(30,33)+linear(39,42)
                if light_subtype=='lambda'and length>=13:return linear(24,27)+linear(29,28+length-8)+linear(39,42)
                tail=5 if length>=15 else 4
                return linear(24,26)+linear(29,28+length-3-tail)+linear(43-tail,42)
            if length<11:
                missing=11-length
                deleted=(set(range(28,28+missing))if scheme=='kabat'else set(range(31,31+missing))if scheme=='chothia'else set(range(31-missing,31)))
                return [str(x)for x in range(24,35)if x not in deleted]
            site=27 if scheme=='kabat' else 30
            return linear(24,site)+[str(site)+string.ascii_uppercase[i]for i in range(length-11)]+linear(site+1,34)
        if region==2:
            if not 7<=length<=15:raise UnsupportedNumbering('L2 supported lengths 7..15')
            if scheme=='imgt':return centered(56,65,length-5)+linear(66,70)
            if scheme=='aho':return linear(58,57+length-4)+linear(69,72)
            site=52 if scheme=='martin'else 54
            return linear(50,site)+[str(site)+string.ascii_uppercase[i]for i in range(length-7)]+linear(site+1,56)
        if region==3:
            if scheme=='imgt':return centered(105,117,length)
            if scheme=='aho':
                if complete and not 2<=length<=32:raise UnsupportedNumbering('AHo L3 exceeds available slots')
                return linear(107,106+length-2)+['137','138'] if complete else linear(107,106+length)
            if not complete:return linear(89,88+length)
            if not 5<=length<=35:raise UnsupportedNumbering('L3 supported lengths 5..35')
            if length<=9:return linear(89,88+length-2)+['96','97']
            return linear(89,95)+['95'+string.ascii_uppercase[i]for i in range(length-9)]+['96','97']
    raise UnsupportedNumbering('Unknown loop')


def framework3_labels(chain,length,scheme):
    """Compute the conserved FR3 path and scheme-specific indel site.

    The positional path depends on segment length, chain, and scheme.
    The C-terminal conserved cysteine stays at its published anchor. IMGT/AHo
    fill their optional central slots before introducing extra insertion labels.
    """
    base=27 if chain=='H'else 32
    if chain=='H'and not 27<=length<=40:raise UnsupportedNumbering('H FR3 exceeds supported path lengths')
    if chain=='L'and not 30<=length<=40:raise UnsupportedNumbering('L FR3 exceeds supported path lengths')
    if scheme in ('imgt','aho'):
        start=(75 if chain=='H'else 71)if scheme=='imgt'else(77 if chain=='H'else 73)
        end=84 if scheme=='imgt'else 86
        tail=104 if scheme=='imgt'else 106
        width=end-start+1+(tail-end)
        if length<=width:return linear(start,end-(width-length))+linear(end+1,tail)
        extra=length-width
        insert=[str(end)+'.'+str(i)for i in range(1,extra+1)]if scheme=='imgt'else[str(end)+string.ascii_uppercase[i]for i in range(extra)]
        return linear(start,end)+insert+linear(end+1,tail)
    start,finish=(66,92)if chain=='H'else(57,88)
    site=(72 if scheme=='martin'else 82)if chain=='H'else(68 if scheme=='martin'else 66)
    count=length-base
    if count<0:raise UnsupportedNumbering('Short light FR3 needs alignment-specific deletion mapping')
    return linear(start,site)+[str(site)+string.ascii_uppercase[i]for i in range(count)]+linear(site+1,finish)


class AnchorNumberer:
    """Number a single variable domain using a fixed set of public frameworks.

    Framework matches and segment boundaries are computed for each input.
    """

    def __init__(self, framework_asset):
        """Load the JSON asset and validate chain, label, and slice consistency."""
        data=json.loads(Path(framework_asset).read_text())
        if data['schema']!='public_numbered_frameworks_v2':raise ValueError('Unsupported public framework schema')
        self.templates=data['templates']
        self.coordinate_paths={}
        for path in data.get('framework_position_paths',{}).get('paths',[]):
            key=(path['chain'],path.get('light_subtype'),path['frame'],tuple(path['labels']['kabat']))
            self.coordinate_paths[key]=path['labels']
        self.framework_profiles=data.get('framework_profiles',{}).get('fragments',[])
        self.profile_groups={}
        self.packed_profiles={profile['id']:pack_residues(profile['sequence'])for profile in self.framework_profiles}
        for profile in self.framework_profiles:
            key=(profile['chain'],profile.get('light_subtype'),profile['frame'])
            self.profile_groups.setdefault(key,[]).append(profile)
        for profile in self.framework_profiles:
            assert set(profile['labels'])==set(SCHEMES)
            assert all(len(labels)==len(profile['sequence'])and all(p.startswith(profile['chain'])for p in labels)for labels in profile['labels'].values())
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
                    refs.append(ordered_framework_alignment(other['sequence'][u:v],query,deletion))
            return max(refs,key=lambda item:item[0])
        if key=='fr4' and len(query)<len(ref):
            # A C-terminal truncation cannot right-align itself to terminal labels.
            return 0.0,[(i,i) for i in range(len(query))]
        return align(ref,query,deletion)if key=='fr3'else ordered_framework_alignment(ref,query,deletion)

    def framework_choice(self, template, key, query):
        """Return a scored mapping from a small public framework-fragment panel.

        Only framework fragments are present in this panel. A donor must match
        chain/locus, cannot leave query residues without labels, and must improve
        the base alignment by the configured evidence margin before replacing
        its positional interpretation.
        """
        a,b=template['segments'][key]
        score,pairs=self.framework_alignment(template,key,query)
        labels={scheme:template['labels'][scheme][a:b]for scheme in SCHEMES}
        base=(score,pairs,labels,'base_template')
        best=base;best_valid=not any(i is None and j is not None for i,j in pairs)
        subtype=template.get('light_subtype','kappa')if template['chain']=='L'else None
        packed_query=pack_residues(query)
        for profile in self.profile_groups.get((template['chain'],subtype,key),[]):
            reference=profile['sequence']
            if len(reference)<len(query):continue
            # FR2 indel hypotheses are searched explicitly by segment length.
            # A longer donor cannot create an insertion label in an ordinary
            # length segment by deleting an unrelated conserved position.
            if key=='fr2'and len(reference)!=len(query):continue
            deletion=[-5.0]*len(reference)
            if key=='fr1'and len(query)<len(reference):
                for i in range(len(reference)-len(query)):deletion[i]=-.1
            if len(reference)==len(query):
                candidate_score=packed_substitution_score(self.packed_profiles[profile['id']],packed_query)
                # Allocate traceback pairs only for a winning equal-length donor.
                candidate_pairs=None
            elif key=='fr4'and len(query)<len(reference):
                candidate_pairs=[(i,i)for i in range(len(query))]
                candidate_score=sum(0 if 'X'in(u,v)else 3 if u==v else -1 for u,v in zip(reference,query))
            else:candidate_score,candidate_pairs=ordered_framework_alignment(reference,query,deletion)
            if candidate_pairs is not None and any(i is None and j is not None for i,j in candidate_pairs):continue
            margin=4 if key in ('fr1','fr4')else 8
            if best_valid and candidate_score<base[0]+margin:continue
            if not best_valid or candidate_score>best[0]:
                if candidate_pairs is None:candidate_pairs=[(i,i)for i in range(len(query))]
                best=(candidate_score,candidate_pairs,profile['labels'],profile['id']);best_valid=True
        return best

    def select_fr4(self, template, sequence, start):
        """Score full FR4-prefix evidence rather than require an exact anchor.

        Conserved Trp/Phe and glycines can mutate. Candidate windows are judged
        against public framework residues; the full-window score favors a real
        framework over a short accidental CDR3 motif. A minimum average score
        rejects weak fragments and leaves truly incomplete domains unmapped.
        """
        a,b=template['segments']['fr4'];reference=template['sequence'][a:b]
        candidates=[]
        for position in range(start,len(sequence)-3):
            query=sequence[position:position+len(reference)]
            score=sum(0 if 'X' in (u,v) else 3 if u==v else -1 for u,v in zip(reference,query))
            if score/len(query)<1.2:continue
            # Very short terminal evidence requires the conserved aromatic start.
            if len(query)<7 and query[0] not in 'FW':continue
            candidates.append(((score,len(query),-(len(sequence)-position-len(query))),position))
        if candidates:return max(candidates)[1]
        # Preserve a recognizable terminal motif when the tail is too divergent
        # for the richer framework-prefix score, rather than lose old coverage.
        pattern=r'W[GAX][A-Z]G'if template['chain']=='H'else r'[FW]G[A-Z]G'
        hits=[start+m.start()for m in re.finditer(pattern,sequence[start:])]
        if not hits:return None
        return max(hits,key=lambda position:sum(0 if 'X'in(a,b)else 3 if a==b else -1 for a,b in zip(reference,sequence[position:position+len(reference)])))

    def locate(self,sequence):
        """Choose a supported chain/template/segmentation hypothesis, or None.

        Enumerate candidate anchors and supported segment lengths, align the
        framework pieces, reject weak matches, then rank the remaining choices.
        Deterministic scores rank the supported segmentation hypotheses.
        Returned c1/c2 are anchor indices; w1 is the candidate start of FR2.
        """
        choices=[];best_score=float('-inf')
        # Candidate FR3 endpoints retain at least two letters of the local YYC
        # motif, allowing Y/F interchange. This hard filter differs from soft C/W bonuses.
        terminals=[i for i in range(2,len(sequence)) if sum((a in 'YF')if b=='Y'else(a=='C') for a,b in zip(sequence[i-2:i+1],'YYC'))>=2]
        for template_id,t in self.templates.items():
            chain=t['chain'];cached={};tail_cache={}
            def framework_score(key,a,b):
                """Memoize repeated slice alignments within this template search."""
                cachekey=(key,a,b)
                if cachekey not in cached:cached[cachekey]=self.framework_choice(t,key,sequence[a:b])[0]
                return cached[cachekey]
            # Search the first 71 offsets for the first anchor. c1 need not be C:
            # the residue contributes a bonus rather than an unconditional veto.
            # c2 anchors FR3; the distances imply candidate loop1/loop2 lengths.
            # Strong conserved-anchor hypotheses are evaluated first. The final
            # tie key still prefers the original increasing c1 order.
            c1_candidates=sorted(range(min(71,len(sequence))),key=lambda i:(sequence[i]!='C',i))
            for c1 in c1_candidates:
                start=max(0,c1-t['c1']);end1=c1+(4 if chain=='H' else 1)
                if end1>len(sequence)or end1<=start:continue
                # Reject weak FR1 evidence once, before exploring downstream lengths.
                fr1_score=framework_score('fr1',start,end1)
                if fr1_score/(end1-start)<.2:continue
                # Every remaining framework residue contributes at most +3;
                # penalties can only decrease the final score. This upper bound
                # safely prunes hypotheses that cannot beat an existing choice.
                remaining_widths=(40,50)if chain=='H'else(45,56)
                maximum_bonus=(6 if sequence[c1]=='C'else 0)+3+(12 if chain=='L'else 0)
                upper=max((fr1_score+3*width+maximum_bonus)/(end1-start+width)for width in remaining_widths)
                if upper<best_score-1e-12:continue
                for loop1length in (range(7,17) if chain=='H' else range(8,21)):
                    w1=end1+loop1length
                    for fr2length in ([13,14] if chain=='H' else [15,16]):
                        afterfr2=w1+fr2length
                        if afterfr2>len(sequence):continue
                        if framework_score('fr2',w1,afterfr2)/fr2length<.2:continue
                        for c2 in terminals:
                            for fr3length in (range(27,37) if chain=='H' else range(30,41)):
                                startfr3=c2+1-fr3length;n2=startfr3-afterfr2
                                if not (15<=n2<=24 if chain=='H' else 7<=n2<=15):continue
                                if c2+1>len(sequence):continue
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
                                    if c2 not in tail_cache:tail_cache[c2]=self.select_fr4(t,sequence,c2+1)
                                    start4=tail_cache[c2]
                                    if start4 is not None:
                                        fr4seq=sequence[start4:start4+10]
                                        # FR4 motifs supply an additional fixed
                                        # +12 kappa/lambda preference when they match;
                                        # contradictory types are rejected; ambiguous tails add none.
                                        lambda_pattern=bool(re.fullmatch(r'.{6}[LIV][TS]V[LV]',fr4seq))
                                        kappa_pattern=bool(re.fullmatch(r'.{7}[ED][ILV][KR]',fr4seq))
                                        is_lambda=t.get('light_subtype')=='lambda'
                                        if lambda_pattern and not is_lambda:continue
                                        if kappa_pattern and is_lambda:continue
                                        if lambda_pattern or kappa_pattern:subtype_support=12
                                # Normalize by framework residues so candidates
                                # of different lengths are compared on this scale.
                                score=(sum(scores)+anchor_support+subtype_support)/sum(lengths)
                                # A shortened FR2 must improve sequence evidence, rather
                                # than win solely by shrinking the normalization divisor.
                                score-=.05*max(0,(14 if chain=='H'else 15)-fr2length)
                                score-=.03*max(0,(30 if chain=='H'else 32)-fr3length)
                                if sum(scores)/sum(lengths)<.45:continue
                                # Broader length support requires stronger joint framework
                                # evidence; calibrated on diagnostic positives/controls.
                                if score<1.3:continue
                                # Prefer the complete framework when support is tied.
                                choices.append((score,chain,segments,c1,c2,w1,template_id))
                                best_score=max(best_score,score)
        if not choices:return None
        return max(choices,key=lambda x:(x[0],sum(b-a for a,b in x[2].values()),x[1],x[6],-x[3]))

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
        # Without a matching FR4, process the remaining tail as incomplete.
        h3start=c2+1
        fr4=self.select_fr4(t,sequence,h3start)
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
                    if key=='fr3'and(chain=='H'or b-a>=32):
                        labels=framework3_labels(chain,b-a,scheme)
                        for j,label in enumerate(labels,a):positions[j]=chain+label
                        continue
                    _,pairs,frame_labels,profile_id=self.framework_choice(t,key,sequence[a:b])
                    # Reconstruct the Kabat framework path after deletions, then
                    # convert that path to each scheme. A deleted Kabat slot
                    # need not have the same physical gap in IMGT or AHo.
                    mapped_kabat=['']*(b-a)
                    for ref_index,query_index in pairs:
                        if ref_index is not None and query_index is not None:mapped_kabat[query_index]=frame_labels['kabat'][ref_index]
                    path_key=(chain,t.get('light_subtype','kappa')if chain=='L'else None,key,tuple(mapped_kabat))
                    converted=self.coordinate_paths.get(path_key)
                    if converted is not None:
                        for j,label in enumerate(converted[scheme],a):positions[j]=label
                        continue
                    for i,j in pairs:
                        if i is None and j is not None:raise UnsupportedNumbering('Framework insertion outside validated anchors')
                        if i is not None and j is not None:positions[a+j]=frame_labels[scheme][i]
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
                if chain=='H' and scheme=='aho' and loops[2][1]-loops[2][0]>=19:
                    positions[loops[2][0]+6:]=['']*(len(sequence)-loops[2][0]-6)
                    notes.append('Legacy abYsis AHo partial output for H2 length>=19; length-conditioned compatibility limit')
            except UnsupportedNumbering as error:
                # Keep the diagnostic reason and input-length correspondence.
                # Do not publish partly assigned labels for a rejected scheme.
                result[scheme]={'status':'unsupported','reason':str(error),'positions':['']*len(sequence)};continue
            result[scheme]={'status':'numbered','positions':positions}
        return {'status':'numbered','chain':chain,'framework_template':template_id,'alignment_score_not_probability':score,'segments':segments,'loops':loops,
                'complete_fr4':complete,'schemes':result,'compatibility_notes':notes}
