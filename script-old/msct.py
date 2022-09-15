class Event:
    
    def __init__ (self, lab, pre, post) :
        pre.sort()
        post.sort()
        self.lab = lab
        self.pre = pre
        self.post = post
        
    def enabled (self, marking) :
        return all([c in marking for c in self.pre]+[not (c in marking) for c in self.post])
        
class ENS:
    
    def __init__ (self, nConditions=0, marking=[], events=[]):
        self.C = [x+1 for x in range(nConditions)]
        marking.sort()
        self.marking = marking
        self.E = events
    
    def new_event (self,label, pre, post) : # modify net by adding a new event
        self.E.append(Event(label, pre, post))
        self.C.extend([c for c in (pre+post) if not (c in self.C) ])
        E.sort(key = lambda e: e.lab)
   
    def set_marking (self, marking) : # set net to given state
        self.C.extend([c for c in marking if not (c in self.C)])
        marking.sort()
        self.marking = marking
    
    def is_EFC(self) : # returns true when the net is extended free choice
        pre_check = all([all([e1.pre==e2.pre for e2 in self.E if any([c in e1.pre for c in e2.pre])]) for e1 in self.E])
        post_check = all([all([e1.post==e2.post for e2 in self.E if any([c in e1.post for c in e2.post])]) for e1 in self.E])
        return (pre_check and post_check)
    
    def enabled(self): # collection of events which are enabled at the current state
        return [e for e in self.E if e.enabled(self.marking)]
    
    def potential(self): # collection of events which are in potential conflict 
        #\{e\in E \mid \exists e'\in E:pre(e)\cap pre(e')\neq \emptyset \lor post(e)\cap post(e')\neq \emptyset\} 
        return [e1 for e1 in self.E if any([any([c in e2.pre for c in e1.pre]+[c in e2.post for c in e1.post]) for e2 in self.E if e1 != e2])]
    
    def conflicts(self): # collection of enabled events which are in a conflict at current state 
        return [e for e in self.enabled() if len(self.conflicting(e))>0]
        
    def conflicting(self, e) : #returns the list of enabled events in conflict with e
        return [e1 for e1 in self.enabled() if e1 != e if any([c in e.pre for c in e1.pre]+[c in e.post for c in e1.post])]
    
    def concurrent(self, e) : # returns the list of concurrently enabled events 
        return [e1 for e1 in self.enabled() if e1 != e if not (any([c in e.pre for c in e1.pre]+[c in e.post for c in e1.post]))]
                          
    def fire (self, e) : #fires the event e
        if e in self.enabled():
            self.marking = [c for c in self.marking if not c in e.pre]+e.post
            self.marking.sort()
        else:
            print('The required event is not enabled.')
            print('Only the following events can fire at current state: '+str([e.lab for e in self.enabled()]))
    
    
    def build_ecuts (self, events) : # recursively builds the maximal sets without enabled conflicts within input the list of events
        steps=[]
        for e1 in events: 
            result = self.build_ecuts([e2 for e2 in events if not (e2 in self.conflicting(e1) or e2==e1)])
            if len(result)!=0: # result contains the steps non conflicting with e1
                for rec in result:
                    new= [e1]+rec
                    new.sort(key = lambda x: x.lab) # to treat list as set
                    if not (new in steps): #avoids repeating equal steps
                        steps.append(new)
            else: # build_ecuts was call with empty list
                steps.append([e1])
        return steps
    
    def get_steps(self): # returns a list of available steps at current marking
        return self.build_ecuts(self.enabled())

class MSCT_node:
    
    def __init__(self, parent, marking, labels=[]):
        self.parent = parent
        self.children = [] #tree strcture
        self.marking = marking
        self.observed = labels
        self.rep_anc = None # ancestor with the same marking
        self.rep_des = [] # descendants with the same marking
        self.footprints = [] # collection of footprints reachable from this node
        self.is_deadlock= False
        self.is_leaf = False
        self.max_runs = []
        self.seen=[]

    def get_children(self, pn):
        s_m = [c for c in self.marking]
        pn.set_marking(s_m)
        steps = pn.get_steps()
        self.children = [] # I couldn'tn figure out why this is needed 
        # self.children behaves like a class attribute?????
        for step in steps:
            lab = [l for l in self.observed]
            for e in step:
                if not( e.lab in lab):
                     lab.append(e.lab)
                pn.fire(e)
            lab.sort()
            child = MSCT_node(self, [c for c in pn.marking], lab)
            self.children.append(child)
            pn.set_marking(s_m)
    
    def build_subtree(self, pn, visited=[]): #returns the list of accessible footprints
        self.get_children(pn)
        if self.children == []: # no children are found iff self is a deadlock
            self.is_deadlock = True
            self.is_leaf = True
            # update footprints
            self.footprints = [[e for e in self.observed]]
            return [[e for e in self.observed]]
        else:
            rep = [anc for anc in visited if anc.marking==self.marking]
            if rep==[]:
                for child in self.children:
                # first check if child is repeated marking
                    fps = child.build_subtree(pn,visited+[self])
                # update footprints
                    self.footprints.extend([fp for fp in fps if not(fp in self.footprints)])
                return self.footprints
            else:
                self.rep_anc=rep[0];
                self.rep_anc.rep_des.append(self)
                self.is_leaf=True
                self.footprints = [[e for e in self.observed]]
                return self.footprints
            
    def delete_all(self): # clean up the mess
        for child in self.children:
            if child.children == []:
                del child
            else:
                delete_all(child)
        del self
        
def build_MSCT(pn): # just for the sake of it
    backup = [c for c in pn.marking]
    root = MSCT_node(None, pn.marking)
    root.footprints = root.build_subtree(pn)
    pn.set_marking(backup)
    return root


def explore(node,k=0, branch=[]): # dfs on the tree
    print('node at depth: '+str(k))
    print('   with marking:'+str(node.marking)) 
    print('   trace:'+str(node.observed))
    print('   on branch: '+str(branch))
    print('   children:'+str([child.marking for child in node.children]), end='')
    if node.is_leaf:
        print('node is a leaf', end='')
        if node.is_deadlock:
            print(' and a deadlock', end='')
    print('')
    b=branch+[node.marking]
    for child in node.children:
        explore(child,k+1,b)
    print('<<<'+str(k))

def reset(node): # dfs on the tree)
    node.max_runs=[]
    for child in node.children:
        reset(child)

def reveals(a, b):
    return all([b in run for run in root.footprints if a in run])

# The following simulates exploration of the full, possibly infinite tree infinite.
# the parameter fp contains the labels of the path explored so far
# A path is considered explored when the same node has been visited twice with the same fp value
def get_footprints(node, fp=[]):
    # get_footprints is called on a node with the current footprint of the path we are exploring
    if node.is_leaf: # If we reach a leaf of the prefix:
        if not(node.is_deadlock) and not(fp in node.seen): # stop condition for recursion
            node.seen.append(fp)
             # node.max_runs is the collection of footprints of all the maximal paths that contain the current node
#            node.max_runs.append(fp) # if fp was not already a maximal run, we add it as so
            for child in node.rep_anc.children:
                fpc=fp+[l for l in child.observed if not(l in fp)]
                fpc.sort()
                rt_fps = get_footprints(child,fpc)
                node.max_runs.extend([run for run in rt_fps if not(run in node.max_runs)])
                # node.rep_anc points to the node which has the same marking as the current one 
        else: 
            if not(fp in node.max_runs):
                node.max_runs.append(fp)            
    else:
        for child in node.children: 
            # If the current node is not a leaf, we simply continue dfs on its children
            fpc=fp+[l for l in child.observed if not(l in fp)]
            fpc.sort()
            rt_fps = get_footprints(child,fpc)
            node.max_runs.extend([run for run in rt_fps if not(run in node.max_runs)])
            # max_runs of current node has all max_runs of all its children
    return node.max_runs

def excludes(a,b):
    return all([not(b in run) for run in root.max_runs if a in run])

#breadth first search, I think this migt be more time-efficient
    
def fpbfs(root):
    fifo = [(root,[])]
    footprints=[]
    while not (fifo==[]):
        (node, fp)= fifo[0]
        fifo=fifo[1:]
        fp=fp+[l for l in node.observed if not(l in fp)] # observed labels are added to fp
        fp.sort() # we sort the labels to avoid spurious repetitions
        if node.is_leaf: # If we reach a leaf of the prefix:
            if not(fp in node.max_runs):
            # node._max runs contains the footprints of the runs which that finally repeat the segment  [node.rep_anc, node]
                if node.is_deadlock:
                    if not(fp in footprints):
                        footprints.append(fp) # if fp was not already a maximal run, we add it as so
                else:
                    fp=fp+[l for l in node.rep_anc.observed if not(l in fp)]
                    fp.sort()
                    if not(fp in footprints):
                        footprints.append(fp) # if fp was not already a maximal run, we add it as so
                    # and continue exploring this path on the ancestor with the same marking
                    node.max_runs.append(fp)
                    fifo.extend([(child, fp) for child in node.rep_anc.children])
        else:
            fifo.extend([(child, fp) for child in node.children])
    return footprints

def fpdfs(root):
    filo = [(root,[])]
    footprints=[]
    while not (filo==[]):
        (node, fp)= filo[0]
        filo=filo[1:]
        fp=fp+[l for l in node.observed if not(l in fp)] # observed labels are added to fp
        fp.sort() # we sort the labels to avoid spurious repetitions
        if node.is_leaf: # If we reach a leaf of the prefix:
            if not(fp in node.max_runs):
            # node._max runs contains the footprints of the runs which that finally repeat the segment  [node.rep_anc, node]
                if node.is_deadlock:
                    if not(fp in footprints):
                        footprints.append(fp) # if fp was not already a maximal run, we add it as so
                else:
                    fp=fp+[l for l in node.rep_anc.observed if not(l in fp)]
                    fp.sort()
                    if not(fp in footprints):
                        footprints.append(fp) # if fp was not already a maximal run, we add it as so
                    # and continue exploring this path on the ancestor with the same marking
                    node.max_runs.append(fp)
                    filo = [(child, fp) for child in node.rep_anc.children]+filo
        else:
            filo=[(child, fp) for child in node.children]+filo
    return footprints


