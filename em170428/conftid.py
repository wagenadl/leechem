import confidence as cf;
def tid_by_confidence(threshold = 50):
    c = cf.Confidence();
    res = [];
    
    for tid in c.sbemconf:
        try:
            if(c.sbemconf[tid]>=threshold):
                res.append(tid);
        except:
            pass;
    return res;