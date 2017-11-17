import pandas as pd
import json
from collections import Counter
from collections import defaultdict
import time
from urllib.request import urlopen
from urllib.request import HTTPError
import requests
import re


def get_class(queries, chunksize = 5):
    in_process = []
    s = requests.Session()
    for i in range(0, len(queries), chunksize):
        to_query = queries[i:i+chunksize]
        label = 'id%d' % i
        payload = {'query_input': '\n'.join(to_query), 'query_type': 'STRUCTURE', 'label': label}
        req = s.post('http://classyfire.wishartlab.com/queries.json', data=json.dumps(payload),
                     headers={'Content-Type': "application/json"})
        if req.status_code != 201:
            print("POST Something went wrong...")
            print(req.status_code)
            print(req.content)
        else:
            in_process.append((label, req.json()['id']))
        #time.sleep(10)
    
    return in_process

def poll(in_process):
    taxon_order = ['kingdom', 'superclass', 'class', 'subclass']
    s = requests.Session()
    records = []
    in_progress_counts = defaultdict(int)
    while in_process:
        id_, query = in_process.pop(0)
        req = s.get("http://classyfire.wishartlab.com/queries/%s.json" % query)
        if req.status_code != 200:
            print("GET something went wrong...")
            continue
        try:
            result = req.json()
        except:
            in_process.append((id_, query))
            in_progress_counts[query] += 1
            continue
        if result['classification_status'] != 'Done':
            if in_progress_counts[query] < 100:
                in_process.append((id_, query))
                in_progress_counts[query] += 1
                continue
        for entry in result['entities']:
            if sum([bool(re.search('report', x)) for x in list(entry.keys())]):
                continue
            if type(entry)==list:
                entry = entry[0]
            if not sum([bool(re.search('identifier', x)) for x in list(entry.keys())]):
                continue
            entry_results = {'#SampleID': entry['identifier']}
            if 'ancestors' in entry:
                for ai, anc in enumerate(entry['ancestors']):
                    anc = anc.replace(' ','_').replace('-','_').replace('/','_').replace(',','_')
                    entry_results[anc] = True
            if 'substituents' in entry:
                for si, subst in enumerate(entry['substituents']):
                    substituent = subst.replace(' ','_').replace('-','_').replace('/','_').replace(',','_')
                    entry_results[substituent] = True
            for taxon in taxon_order:
                if taxon in entry:
                    if entry[taxon] is not None:
                        entry_results['%s_id' % taxon] = entry[taxon].get('chemont_id', None)
                        entry_results['%s_name' % taxon] = entry[taxon].get('name', None)
            #for ci, ch in enumerate(entry['predicted_chebi_terms']):
            #    entry_results['predicted_chebi_term_%d' % ci] = ch
            if 'molecular_framework' in entry:
                entry_results['molecular_framework'] = entry['molecular_framework']
            records.append(entry_results)
        #time.sleep(10)
    return pd.DataFrame(records) 

def query_inchikey(keylist):
    count_failed = 0
    dlist = []
    dmetadatalist = []
    for i in range(len(keylist)):
        key = keylist[i]
        url = 'http://classyfire.wishartlab.com/entities/' + key + '.json' 
    
        try:
            r = urlopen(url) 
            fn = r.read().decode('utf8')
        except  HTTPError as err:
            fn = ''
        #try:
        if fn != '':
            mdict = {}
            data = json.loads(fn) 

            mdict['inchikey'] = key 
            if sum([bool(re.match('kingdom', x)) for x in data.keys()]) > 0 and data['class'] is not None:
                mdict['kingdom'] = data['kingdom']['name']
            if sum([bool(re.match('superclass', x)) for x in data.keys()]) > 0 and data['superclass'] is not None :
                mdict['superclass'] = data['superclass']['name']
            if sum([bool(re.match('class', x)) for x in data.keys()]) > 0 and data['class'] is not None :
                mdict['class'] = data['class']['name']
            if sum([bool(re.match('subclass', x)) for x in data.keys()]) > 0 and data['subclass'] is not None:
                mdict['subclass'] = data['subclass']['name']
            if sum([bool(re.match('direct_parent', x)) for x in data.keys()]) > 0 and data['direct_parent'] is not None:
                mdict['direct_parent'] = data['direct_parent']['name']
            if sum([bool(re.match('molecular_framework', x)) for x in data.keys()]) > 0 and data['molecular_framework'] is not None:
                mdict['molecular_framework'] = data['molecular_framework']
            dmetadatalist.append(mdict)
        #except HTTPError as err:
        else:
            count_failed += 1
            dmetadatalist.append({})
    
    df_metares = pd.DataFrame.from_dict(dmetadatalist)  
    return(df_metares)
