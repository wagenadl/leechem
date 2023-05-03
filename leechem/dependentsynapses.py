#!/usr/bin/python3

import sbemdb
import re

r = re.compile(':(\d+)')

db = sbemdb.SBEMDB()

nsoma = db.nodeDetails('tid==444 and typ==1')
for a in nsoma:
    somanid = a

c = db.db.cursor()
def recurse(nid, seen={}, tags=[]):
    seen = seen.copy()
    tags = tags.copy()
    seen[nid] = 1
    c.execute('select typ from nodes where nid==%i' % nid)
    typ = c.fetchone()[0]
    if typ==6:
        if len(tags):
            print('postsyn', nid, tags)
            c.execute('''select tag from tags where nid==%i''' % nid)
            newtag = "; ".join(tags)
            tgs = c.fetchall()
            if len(tgs):
                oldtag = tgs[0][0]
                tag = oldtag + "; " + newtag
                c.execute('update tags set tag="%s" where nid==%i' % (tag, nid))
            else:
                tag = newtag
                c.execute('insert into tags (nid, tag) values (%i,"%s")'
                          % (nid, tag))
    else:
        c.execute('''select tag,x,y,z from tags inner join nodes
        on tags.nid==nodes.nid where nodes.nid==%i''' % nid)
        tgs = c.fetchall()
        for tg in tgs:
            if 'uncertain' in tg[0].lower():
                x = tg[1] *.0055
                y = tg[2] *.0055
                z = tg[3] *.050
                tag = tg[0]
                m = r.search(tag)
                if m:
                    tag = 'b:' + m.group(1)
                tags.append(tag + (' from %.2f,%.2f,%.2f' % (x,y,z)))
    c.execute('select nid2 from nodecons where nid1==%i' % nid)
    nnn = c.fetchall()
    for nn in nnn:
        n = nn[0]
        if n not in seen:
            recurse(n, seen, tags)
            
recurse(somanid)    
db.db.commit()
