#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, re, sys
from pathlib import Path

src=Path(sys.argv[1]); data=src.read_bytes(); text=data.decode('utf-8',errors='replace')
expected=[5,6,8,52,63,64,65,66,67,69,71,73,74,75,79,80,81,82,88,89,90,93,94,98]
symbols=sorted(set(re.findall(r'(?<![A-Za-z0-9_$`])([A-Za-z$][A-Za-z0-9$`]*)\s*\[',text)))
selected=[]
for n,line in enumerate(text.splitlines(),1):
    if any(k in line.lower() for k in ('diagram','unsolved','mass relation','zwmatrix','zw matrix')):
        selected.append(f'{n}: {line[:3000]}')
mentions=sorted(set(int(x) for x in re.findall(r'(?i)(?:diagram|case)\s*(?:no\.?\s*)?(\d{1,3})',text) if 1<=int(x)<=200))
out={
 'sha256':hashlib.sha256(data).hexdigest(), 'bytes':len(data),
 'prefix':text[:120], 'Notebook_count':text.count('Notebook['),
 'Cell_count':text.count('Cell['), 'BoxData_count':text.count('BoxData['),
 'symbol_count':len(symbols),'symbols':symbols,
 'diagram_mentions':mentions,'expected_unresolved_24':expected,
 'expected_mentions_present':[x for x in expected if x in mentions],
 'selected_line_count':len(selected),
 'scope':'Source inventory only; no mathematical conclusion without exact execution or an independent port.'}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
Path(sys.argv[3]).write_text('\n'.join(selected)+'\n')
