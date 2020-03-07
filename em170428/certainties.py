from sbemdb import SBEMDB, clean_db;

def decode_tags(tags):
    tags_array = tags.split(';');
    if(len(tags_array) == 1):
        tags_array = tags.split('.');
    if(len(tags_array) == 1):
        tags_array = tags.split(',');
    tags_dictionary = {"s":100,"b":100};
    for tag in tags_array:
        tag_key_value = tag.split(':');
        if(len(tag_key_value) == 2):
            try:
                tags_dictionary[tag_key_value[0].strip().lower()] = int(tag_key_value[1]);
            except:
                try:
                    tags_dictionary[tag_key_value[0].strip().lower()] = int(tag_key_value[1].split(" ")[0]);
                except:
                    tags_dictionary[tag_key_value[0].strip().lower()] = tag_key_value[1];
    return tags_dictionary;

def get_tags(nid, cursor):
    query = cursor.execute(f'select tag from tags where nid ={nid}');
    result = query.fetchall();
    if(len(result) != 0):
        tags = decode_tags(result[0][0]);
    else:
        tags = decode_tags("");
    
    b = tags['b'];
    s = tags['s'];
    return b,s;

def uncertainties():
    db = SBEMDB()
    db = clean_db(db)
    cursor = db.db.cursor()
    (xx, yy, zz, pretid, posttid, synid, prenid, postnid) = db.synapses(extended=True);
    uncertain = dict();
    for i in range(0,len(synid)):
        prenidr = prenid[i];
        postnidr = postnid[i];
        prenid_b, prenid_s = get_tags(prenidr, cursor);
        postnid_b, postnid_s = get_tags(postnidr, cursor);
        
        s = prenid_s;
    
        if(prenid_s == 100 or postnid_s == 100):
            s = prenid_s + postnid_s - 100;
            
        if(prenid_s > 1 and postnid_s > 1 and prenid_s != postnid_s):
            uncertain[synid[i]] = 1;
        else:
            uncertain[synid[i]] = prenid_b * postnid_b * s / 100**3;
            
    return uncertain;