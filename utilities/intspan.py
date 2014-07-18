"""
This module provides methods for common operations on integer spans

Author: Readman Chiu rchiu@bcgsc.ca
"""
def cardinality(span):
    """Returns cardinality (size) of span"""
    return abs(int(span[1])-int(span[0])+1)

def overlap(span1, span2):
    """Returns size of overlap of 2 spans
    False if 0
    """
    if int(span1[0]) <= int(span2[1]) and int(span1[1]) >= int(span2[0]):
        if subsume(span1,span2):
            return cardinality(span1)
        elif subsume(span2,span1):
            return cardinality(span2)
        elif int(span1[0]) > int(span2[0]):
            return int(span2[1]) - int(span1[0]) + 1
        else:
            return int(span1[1]) - int(span2[0]) + 1
    else:
        return False

def subsume(span1, span2):
    """Determines if span1 is subsumed in span2"""
    if int(span1[0]) >= int(span2[0]) and int(span1[1]) <= int(span2[1]):
        return True
    else:
        return False

def cardinality_multi(spans):
    """Returns size of overlap of multiple spans"""
    sum_span = 0
    for span in spans:
        sum_span += int(span[1]) - int(span[0]) + 1
        
    return sum_span

def same_blocks(span1, span2):
    """Determines if 2 alignment blocks are the 'same'"""
    a_start = -1
    b_start = -1
    #find first matching block
    for i in range(len(span1)):
        for j in range(len(span2)):
            if subsume(span1[i],span2[j]) or subsume(span2[j],span1[i]):
                a_start = i
                b_start = j
                break

        if a_start >= 0 or b_start >= 0:
            break
        
    if a_start < 0 or b_start < 0:
        return False
    elif a_start > 0 and b_start > 0:
        return False
    elif b_start == len(span2)-1 or a_start == len(span1)-1:
        return False
    else:
        a_end = a_start
        b_end = b_start

        if len(span1) < len(span2):
            j = b_start + 1
            for i in range(a_start+1, len(span1)):
                if subsume(span1[i],span2[j]) or subsume(span2[j],span1[i]):
                    a_end = i
                    b_end = j
                    j += 1
                else:
                    break
                if j == len(span2):
                    break                      
        else:
            i = a_start + 1
            for j in range(b_start+1, len(span2)):
                if subsume(span1[i],span2[j]) or subsume(span2[j],span1[i]):
                    a_end = i
                    b_end = j
                    i += 1
                else:
                    break
                if i == len(span1):
                    break

        if (a_end == len(span1)-1 and a_start == 0) or (b_end == len(span2)-1 and b_start == 0):
            return True
        else:
            return False

def intersect(span1, span2):
    """Intersects 2 list of spans(e.g. exons and alignment blocks) and returns size of overlap"""
    olap_total = 0   
    for span11 in span1:
        for span22 in span2:
            olap = overlap(span11, span22)

            if olap:
                olap_total += olap

    if olap_total > 0:
        return olap_total
    else:
        return False

def union(spans):
    """Unions give list of spans"""
    span_union = []
    for span in spans:
        for span1 in span:
            span1[0] = int(span1[0])
            span1[1] = int(span1[1])
    span_union = spans[0]

    for i in range(1,len(spans)):      
        for bi in range(len(spans[i])):
            olaps = []
            olap_bounds = []

            for bu in range(len(span_union)):             
                if overlap(spans[i][bi], span_union[bu]):
                    olaps.append(bu)
                    olap_bounds.extend(span_union[bu])

            if olaps:
                olap_bounds.extend(spans[i][bi])
                new_span = [min(olap_bounds), max(olap_bounds)]

                olaps.reverse()
                for idx in olaps:
                    span_union.pop(idx)

                span_union.append(new_span)
                
                del olaps[:]
                del olap_bounds[:]
            else:
                span_union.append(spans[i][bi])
            
            span_union.sort(lambda x,y: x[0]-y[0])

    return span_union

def subtract(span1, spans2):
    """Subtracts spans2 from span1
    BEWARE: will modify spans2
    """
    span2_union = union(spans2)

    for span2 in span2_union:
        if subsume(span1, span2):
            return False

    if int(span1[0]) < int(span2_union[0][0]):
        span2_union.insert(0, [int(span1[0])-2, int(span1[0])-1])
    if int(span1[1]) > int(span2_union[-1][1]):
        span2_union.append([int(span1[1])+1, int(span1[1])+2])

    span2_reverse = []
    for i in range(len(span2_union)-1):
        s1 = span2_union[i]
        s2 = span2_union[i+1]
        span2_reverse.append([s1[1]+1, s2[0]-1])

    return span2_reverse

def merge_blocks(blocks):
    """Merges blocks if end coord of first block and start coord of second block are continuous
    Changes input blocks
    """
    merge = []
    start = 0
    while start < len(blocks)-1:
        same = [start]
        for j in range(start+1, len(blocks)):
            if blocks[j][0] - blocks[same[-1]][1] == 1:
                same.append(j)
            else:
                break
        start = same[-1] + 1
        if len(same) > 1:
            merge.append(same)

    merge.reverse()
    for i in merge:
        blocks[i[0]:i[-1]+1] = [[blocks[i[0]][0],blocks[i[-1]][1]]]

