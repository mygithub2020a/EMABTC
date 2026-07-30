#!/usr/bin/env python3
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json,re,sys

@dataclass
class Node:
    kind:str
    value:object
    args:list|None=None

class Parser:
    def __init__(self,s): self.s=s; self.i=0
    def ws(self):
        while self.i<len(self.s) and self.s[self.i].isspace(): self.i+=1
    def parse(self):
        self.ws()
        if self.i>=len(self.s): raise ValueError('eof')
        c=self.s[self.i]
        if c=='"': return Node('str',self.string())
        if c=='{':
            self.i+=1; xs=[]; self.ws()
            if self.s[self.i:self.i+1]=='}': self.i+=1; return Node('list',None,xs)
            while True:
                xs.append(self.parse()); self.ws()
                if self.s[self.i:self.i+1]==',': self.i+=1; continue
                if self.s[self.i:self.i+1]=='}': self.i+=1; break
                raise ValueError(('list',self.i,self.s[self.i:self.i+50]))
            return Node('list',None,xs)
        j=self.i
        while self.i<len(self.s) and self.s[self.i] not in '[],{}' and not self.s[self.i].isspace(): self.i+=1
        atom=self.s[j:self.i]; self.ws()
        if self.s[self.i:self.i+1]=='[':
            self.i+=1; args=[]; self.ws()
            if self.s[self.i:self.i+1]==']': self.i+=1; return Node('func',atom,args)
            while True:
                args.append(self.parse()); self.ws()
                if self.s[self.i:self.i+1]==',': self.i+=1; continue
                if self.s[self.i:self.i+1]==']': self.i+=1; break
                raise ValueError(('func',atom,self.i,self.s[self.i:self.i+80]))
            return Node('func',atom,args)
        return Node('atom',atom)
    def string(self):
        self.i+=1; out=[]
        while self.i<len(self.s):
            c=self.s[self.i]; self.i+=1
            if c=='"': return ''.join(out)
            if c=='\\':
                d=self.s[self.i]; self.i+=1
                if d=='\n': continue
                if d=='\r':
                    if self.i<len(self.s) and self.s[self.i]=='\n': self.i+=1
                    continue
                out.append({'n':'\n','r':'\r','t':'\t','"':'"','\\':'\\'}.get(d,'\\'+d))
            else: out.append(c)
        raise ValueError('unterminated string')

SPECIAL={'\\[IndentingNewLine]':'\n','\\[NewLine]':'\n','\\[InvisibleSpace]':'','\\[Rule]':'->','\\[RuleDelayed]':':>','\\[Equal]':'==','\\[NotEqual]':'!=','\\[LessEqual]':'<=','\\[GreaterEqual]':'>=','\\[And]':'&&','\\[Or]':'||','\\[Function]':'&','\\[Times]':'*','\\[CenterDot]':'*','\\<':'','\\>':''}
def clean(x):
    for a,b in SPECIAL.items(): x=x.replace(a,b)
    return x

def render(n):
    if n.kind in ('str','atom'): return clean(str(n.value))
    if n.kind=='list': return '{'+','.join(render(x) for x in n.args)+'}'
    h=n.value; a=n.args
    if h=='RowBox': return ''.join(render(x) for x in (a[0].args if a and a[0].kind=='list' else a))
    if h in ('BoxData','StyleBox','TagBox','FormBox','TooltipBox','PaneBox','FrameBox'): return render(a[0]) if a else ''
    if h=='InterpretationBox': return render(a[1]) if len(a)>1 else render(a[0])
    if h=='SuperscriptBox': return '('+render(a[0])+')^('+render(a[1])+')'
    if h=='SubscriptBox': return 'Subscript['+render(a[0])+','+render(a[1])+']'
    if h=='FractionBox': return '('+render(a[0])+')/('+render(a[1])+')'
    if h=='SqrtBox': return 'Sqrt['+render(a[0])+']'
    if h in ('GridBox','Column','TemplateBox'): return render(a[0])
    return h+'['+','.join(render(x) for x in a)+']'

def match(s,start):
    depth=0; i=start; ins=False; esc=False
    while i<len(s):
        c=s[i]
        if ins:
            if esc: esc=False
            elif c=='\\': esc=True
            elif c=='"': ins=False
        else:
            if c=='"': ins=True
            elif c=='[': depth+=1
            elif c==']':
                depth-=1
                if depth==0: return i
        i+=1
    raise ValueError('unmatched bracket')

def main():
    if len(sys.argv)!=3: raise SystemExit('usage: extract_input_cells.py INPUT.nb OUTDIR')
    text=Path(sys.argv[1]).read_text(errors='replace'); cells=[]; pos=0
    while True:
        j=text.find('Cell[',pos)
        if j<0: break
        a=j+len('Cell'); b=match(text,a); block=text[j:b+1]
        if re.search(r'\],\s*"Input"(?:,|\])',block,re.S) and 'BoxData[' in block:
            k=block.find('BoxData[')+len('BoxData'); e=match(block,k)
            cells.append(render(Parser(block[k+1:e]).parse()))
        pos=b+1
    out=Path(sys.argv[2]); out.mkdir(parents=True,exist_ok=True)
    for i,c in enumerate(cells,1): (out/f'input_{i:02d}.wl').write_text(c+'\n')
    (out/'manifest.json').write_text(json.dumps({'count':len(cells),'lengths':[len(c) for c in cells]},indent=2)+'\n')
    if len(cells)<4 or 'FindOrder' not in cells[2]: raise SystemExit('Algorithm-II extraction canary failed')
    print(json.dumps({'count':len(cells),'lengths':[len(c) for c in cells]}))
if __name__=='__main__': main()
