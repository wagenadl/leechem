#!/usr/bin/python3

import numpy as np
from . import sbemdb

class Segment:
    '''SEGMENT - Representation of a segment of a tree

    Some terminology:

    - A NODE is a node of the neurite tree with an ID that matches the sbemdb.
    - An EDGE is a connection between adjacent nodes, corresponding to a 
      "nodecon" in the sbemdb, but without its own ID.
    - The FANOUT of a NODE is the number of EDGES that meet there.
    - A SEGMENT is a linear chain of NODES connected by EDGES such that
      all the internal NODES have FANOUT equal to two, whereas the first
      and last NODE have FANOUT not equal to two.

    Segments are organized in a TREE. The ROOT of the tree typically is
      the segment that contains the soma. 
    - CHILD segments by definition connect to the PARENT segment with their
      first node attached to the parent's last node.
    - A TERMINAL SEGMENT is a SEGMENT without children.

    Segments have the following properties:
    - SID - An (almost) arbitrary ID for this segment. (By construction, 
      the root of the tree has ID 0, and a child segment always has a higher
      numbered ID than its parent.)
    - PARENT - The ID of the parent of the segment, or None for the root
    - CHILDREN - A list of the IDs of the children of this segment
    - PATHLEN - The length of the segment, measured as the sum of the 
      lengths of the edges between its nodes
    - HASSYN - True or False depending on whether the segment contains
      at least one synapse (including at its initial or terminal node).
    - DEPTH - Number of segments to traverse to get to the root segment.
      (In future, we'll also measure TRUNKDEPTH, the number of segments
      to traverse to get to a segment that is part of the trunk.)
    - NEIGHBOR_DIST - Euclidean distance between the terminal node of
      this segment and either the final segment of its parent, or the
      first segment of one of its siblings.
    - SKELETON_DIST - A tuple containing: (1) Euclidean distance between 
      the terminal node of this segment and the nearest point on the
      tree that is not on this segment and (2) the "scid" of the edge that
      contains that point (see the "nodecons" table in the sbemdb).
      (SKELETON_DIST is only calculated for terminal segments.)
    '''
    def __init__(self, sid=None, parent=None, nodes=[]):
        '''Construct an empty segment'''
        self.sid = sid
        self.parent = parent
        self.nodes = nodes
        self.children = []
        self.pathlen = None
        self.hassyn = None
        self.depth = 0 # Number of junctions from root
        self.neighbor_dist = None 
        self.skeleton_dist = None
        
    def is_terminal(self):
        '''IS_TERMINAL - Returns true if the segment has no children'''
        return len(self.children)==0

    def edge_count(self):
        '''EDGE_COUNT() returns the number of EDGES in the chain'''
        N = len(self.nodes)
        if N==0:
            return 0
        else:
            return N-1
    
    def __repr__(self):
        r = f'Segment({self.sid})'
        if len(self.nodes)==1:
            r += f' with 1 node: [{self.nodes[0]}]'
        else:
            r += f' with {len(self.nodes)} nodes'
            if len(self.nodes)>=2:
                r += f': [{self.nodes[0]} ... {self.nodes[-1]}]'
        if len(self.children)==1:
            r += f' and 1 child'
        elif len(self.children)==0:
            r += f' without children'
        else:
            r += f' and {len(self.children)} children'
        return r

class Tree(dict):
    '''TREE - Representation of a neuritic tree in terms of segments
    
    A TREE is a dict mapping "segment IDs" to segments.
    The root of the tree has ID zero. 

    - DROP: Remove a terminal segment from the tree and remodel as needed.
'''
    
    def __init__(self, rootnid, db):
        '''Tree(rootnid, db) construct a tree starting from the given node.
        DB must be an SBEMDB'''
        self.seen = set()
        tid = db.fetch(f'select tid from nodes where nid=={rootnid}')[0][0]
        self.tid = tid
        self.rootnid = rootnid
        rows = db.fetch(f'''select nid1, nid2 from nodecons 
        inner join nodes on nid2==nid
        where tid=={tid}''')
        self.cons = {}
        for row in rows:
            if row[0] not in self.cons:
                self.cons[row[0]] = []
            self.cons[row[0]].append(row[1])
        self._explore([rootnid], None)
        del self.cons
        del self.seen
        self._typify(db, tid)
        self._establish_depths()
        self._establish_lengths(db, rootnid)
        self._establish_neighbor_dists(db)
        self._establish_skeleton_dists(db)

    def _idgen(self):
        return len(self)

    def _explore(self, chain, parent):
        nid = chain[-1]
        for n in chain:
            self.seen.add(n)
        while True:
            nids = [n for n in self.cons[nid] if n not in self.seen]
            for n in nids:
                self.seen.add(n)
            if len(nids)==0:
               myid = self._idgen()
               s = Segment(myid, parent, chain)
               self[myid] = s
               return myid
            elif len(nids)==1:
                chain.append(nids[0])
                nid = nids[0]
            else:
               myid = self._idgen()
               s = Segment(myid, parent, chain)
               self[myid] = s
               for n in nids:
                   self[myid].children.append(self._explore([nid, n], myid))
               return myid

    def _establish_depths(self):
        sids = list(self.keys())
        sids.sort()
        for sid in sids:
            pid = self[sid].parent
            if pid is None:
                self[sid].depth = 0
            else:
                par = self[pid]
                self[sid].depth = self[pid].depth + 1

    def _establish_lengths(self, db, rootnid):
        dd = db.distanceAlongTree(rootnid)
        for sid, seg in self.items():
            self[sid].pathlen = dd[seg.nodes[-1]] - dd[seg.nodes[0]]
                   
    def _typify(self, db, tid):
        rows = db.fetch(f'select nid, typ from nodes where tid=={tid}')
        typs = { row[0]: row[1] for row in rows }
        for sid, seg in self.items():
            hassyn = False
            for n in seg.nodes:
                if typs[n]==6:
                    hassyn = True
                    break
            self[sid].hassyn = hassyn

    def is_terminal(self, sid):
        '''IS_TERMINAL(id) returns true if the segment has no children'''
        return len(self[sid].children)==0

    def siblings(self, sid):
        '''SIBLINGS(sid) returns a list of siblings of a given segment.'''
        pid = self[sid].parent
        if pid is None:
            return []
        return [ s for s in self[pid].children if s != sid ]
            
    def drop(self, sid):
        '''DROP(id) removes the named segment from the tree.
        Only terminal segments can be removed.
        The tree is remodeled such that if the removed segment had 
        precisely one sibling, that sibling is concatenated to its 
        parent and removed from the tree.'''
        if not self.is_terminal(sid):
            raise ValueError(f'Cannot drop internal segment {sid}')
        pid = self[sid].parent
        sibs = self.siblings(sid)
        self[pid].children.remove(sid)
        del self[sid]

        def reestablish_depth(tree, sid, d):
            for s in tree[sid].children:
                tree[s].depth = d + 1
                reestablish_depth(tree, s, d + 1)

        if len(sibs)==1:
            sib = sibs[0]
            self[pid].nodes += self[sib].nodes[1:]
            if self[sib].hassyn:
                self[pid].hassyn = True
            self[pid].children = self[sib].children
            self[pid].pathlen += self[sib].pathlen
            for c in self[pid].children:
                self[c].parent = pid
            del self[sib]
            reestablish_depth(self, pid, self[pid].depth)

    def _all_edges(self, db, tid=None):
        '''_ALL_EDGES - Return all the edges in a tree. They are returned
        as a dict from nodecon ID to a pair of node IDs.'''
        if tid is None:
            tid = self.tid
        rows = db.fetch(f'''select ncid, nid1, nid2 from nodecons 
        left join nodes on nodecons.nid1==nodes.nid
        where tid=={tid} and nid1<nid2''')
        edges = {}
        for row in rows:
            edges[row[0]] = ( row[1], row[2] )
        return edges

    def _all_points(self, db, tid=None):
        if tid is None:
            tid = self.tid
        (xx, yy, zz, nids) = db.nodexyz(f'tid=={tid}')
        K = len(xx)
        nodes = {}
        for k in range(K):
            nodes[nids[k]] = ( xx[k], yy[k], zz[k] )
        return nodes

    def _skeleton_dist(self, sid, db, edg=None, pts=None):
        if edg is None:
            edg = self._all_edges(db)
        if pts is None:
            pts = self._all_points(db)
        ncid0 = None
        dist0 = 1e20
        nn = self[sid].nodes
        nid_term = nn[-1]
        p = pts[nid_term]
        for ncid, nids in edg.items():
            n1, n2 = nids
            if n1 not in nn or n2 not in nn:
                # This edge is not part of the named segment
                p1 = pts[n1]
                p2 = pts[n2]
                d1 = sbemdb.LineSegmentGeom.pointDistance(p, p1, p2)
                if d1<dist0:
                    dist0 = d1
                    ncid0 = ncid
        return (dist0, ncid0)
    
    def _establish_skeleton_dists(self, db):
        '''The skeleton distance is the shortest distance between the
        terminal node of a segment and any other point in the tree.'''
        edg = self._all_edges(db)
        pts = self._all_points(db)
        print('Establishing skeleton distances...')
        for s in self.keys():
            if self[s].is_terminal():
                self[s].skeleton_dist = self._skeleton_dist(s, db, edg, pts)
                print(s, end='\r')
            else:
                self[s].skeleton_dist = None
        print()
            
    def _neighbor_dist(self, sid, db, pts=None):
        '''NEIGHBOR_DIST - Distance between terminal node and neighboring edges
        NEIGHBOR_DIST(seg_id, db) returns the distance (in microns) between
        the terminal node of the given chain and the edges in the chain's
        parent or siblings that are adjacent to the given chains base
        node. Second argument must be a SBEMDB to use for measurements.'''
        nid_term = self[sid].nodes[-1]
        nid_base = self[sid].nodes[0]
        if pts is not None:
            p = pts[nid_term]
            p1 = pts[nid_base]
        par_id = self[sid].parent
        if par_id is None:
            return None
        sib_ids = self.siblings(sid)
        nid_neighbors = [ self[par_id].nodes[-2] ]
        for s in sib_ids:
            nid_neighbors.append(self[s].nodes[1])
            if pts is None:
                dd = [ db.nodeToEdgeDistance(nid_term, nid_base, n_nei)
                       for n_nei in nid_neighbors ]
            else:
                dd = [ sbemdb.LineSegmentGeom.pointDistance(p, p1, pts[n_nei])
                       for n_nei in nid_neighbors ]
                
        return np.min(dd)

    def _establish_neighbor_dists(self, db):
        pts = self._all_points(db)
        print('Establishing neighbor distances...')
        for s in self.keys():
            self[s].neighbor_dist = self._neighbor_dist(s, db, pts)

    def prune_stubs(self, db, max_dist_um, max_edge_cnt=1, use_skel=True):
        '''PRUNE_STUBS - Remove very short terminal segments and remodel
        PRUNE_STUBS(db, max_dist_um) prunes terminal segments that have
        SKELETON_DIST less than MAX_DIST_UM and EDGE_COUNT equal to one.
        Optional argument MAX_EDGE_CNT overrides the latter threshold.
        Optional argument USE_SKEL can be set to FALSE to use NEIGHBOR_DIST
        instead of SKELETON_DIST.'''
        cands = set()
        if use_skel:
            def is_short(seg):
                return seg.skeleton_dist[0] <= max_dist_um
        else:
            def is_short(seg):
                return seg.neighbor_dist <= max_dist_um
        for sid, seg in self.items():
            if (seg.is_terminal()
                and seg.edge_count() <= max_edge_cnt
                and is_short(seg)):
                cands.add(sid)
                # CANDS are terminal segments we are willing to prune
                # Careful: If _all_ siblings at some junction are prunable,
                # we should keep the longest one
        keep = set()
        for sid in cands:
            sibs = self.siblings(sid)
            alldel = True
            for s in sibs:
                if s not in cands:
                    alldel = False
                    break
            if alldel:
                sibs.append(sid)
                lens = [ self[s].pathlen for s in sibs ]
                keep.add(sibs[np.argmax(lens)])
        for s in keep:
            cands.remove(s)
            cands = list(cands)
            cands.sort()
        for s in cands[-1::-1]:
            if s in self:
                self.drop(s)
                if not use_skel:
                    self._establish_neighbor_dists(db)
        self._establish_neighbor_dists(db)
        self._establish_skeleton_dists(db)

    def plot(self, db):
        import matplotlib.pyplot as plt
        def xyz(nodes):
            xx=[]
            yy=[]
            zz=[]
            for n in nodes:
                x,y,z = db.onenodexyz(n)
                xx.append(x)
                yy.append(y)
                zz.append(z)
            return xx,yy,zz
        plt.figure()
        for sid, seg in self.items():
            xx,yy,zz = xyz(seg.nodes)
            if seg.depth==0:
                c = 'r'
            elif seg.depth==1 and not seg.is_terminal():
                c = 'b'
            elif seg.is_terminal():
                c = 'k'
            else:
                c = 'y'
                plt.plot(xx,zz, c)
                

    def save_csv(self, ofn):
        with open(ofn, 'w') as fd:
            fd.write('segment_id,branch_id,type,depth,is_synapse,'
                     + 'node_id,point_node_id\n')
            for sid, seg in self.items():
                if seg.depth==0:
                    typ = 'main'
                elif seg.depth==1 and not seg.is_terminal():
                    typ = 'root'
                elif seg.is_terminal():
                    typ = 'ts'
                else:
                    typ = 'is'
                for k in range(1, len(seg.nodes)):
                    fd.write(f'{sid},0,{typ},{seg.depth},{seg.hassyn},'
                             + f'{seg.nodes[k]},{seg.nodes[k-1]}\n')
                    
