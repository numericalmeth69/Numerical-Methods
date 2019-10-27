/* Standard blocker, part of Pogo - a finite element package to 
simulate elastic wave propagation on the GPU
Copyright (C) 2013 Peter Huthwaite

If you find Pogo useful in your academic work, please cite the relevant papers;
information on our latest papers is available at <http://www.pogo-fea.com/>.

This file is part of the blocker, and part of Pogo.

Pogo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Pogo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Pogo.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <list>

#include "block.h"
#include "node.h"
#include <cmath>

#include "../common/err.h"
#include "../common/params.h"

//this is basically an extension of block.cpp
//does the memory mapping within each block
//and calculates the link array to point between blocks

using namespace std;


int block::memArrange(int * nodeBounds, vector<block> &blocks) {
    //return -10 means we can improve things by splitting
    //return -2 means there are problems with the linked blocks (i.e. not ready yet)

    if (blockFlag >= 2) {
        // don't do anything because memory for block is already done
        return 0;
    }

    if (size() == 0) {
        stringstream errSs;
        errSs << "Block has zero size. Block " << blockNum << endl;
        err newError("Blocker - memory arrangement", 2409, errSs.str());
        throw newError;
    }

    if (boundNodes.size()  > MAXBOUNDSIZE) {
        cout << "Block has too many boundary nodes: " << boundNodes.size() << endl;
        return -10;
    }

    if (blockFlag < 1) {
        stringstream errSs;
        errSs << "blockFlag < 1; block should have final shape. Block " << blockNum << endl;
        err newError("Blocker - memory arrangement", 2410, errSs.str());
        throw newError;
    }

    //get a list of the adjacent blocks
    list <int> attachedBlocks;
    getLinks(attachedBlocks);
    //check that the shapes of the attached blocks is finalised
    for (list<int>::iterator p = attachedBlocks.begin();
         p != attachedBlocks.end(); p++) {
        //cout << "Attached: " << *p << " of flag " << blocks[*p].blockFlag << endl;
        if (blocks[*p].blockFlag < 1) {
            //cout << "Error - linked block " << *p << " has not had its shape finalised." << endl;
            return -2;
        }
    }

    //cout << "Getting starting node." << endl;
    //1. get a starting node
    int startNode = -1;
    list <int> bothBounds; //bothBounds is used when a block lies on the boundary (internal or external) of the domain
    for (list<node>::iterator p = boundNodes.begin();
         p != boundNodes.end(); p++) { // loop through the boundary of the block
        int nLinks = 0;
        for (int linkCnt = 0; linkCnt < nNodesLinked[p->nodeNum]; linkCnt++) {
            //loop through the links counting up the number of boundary linked nodes
            //recording one of the attached blocks
            int linkNode = linkMat[(p->nodeNum)*nMaxLinkedNodes+linkCnt];
            if (nodeBlocks[linkNode] == blockNum && flags[linkNode] >= 2) {
                nLinks++;
            }
        }
        if (nLinks == 2 && startNode == -1) {
            startNode = p->nodeNum;
        }
        if (nLinks == 1 && nodeBounds[p->nodeNum] != 0) {
            //this is only linked to one other block boundary
            //and it is on the boundary of the mesh
            //therefore it would form the 'end' of a loop going around the block boundary
            bothBounds.push_back(p->nodeNum);
        }
    }


    if (bothBounds.size() > 0) {
        startNode = bothBounds.front();
        bothBounds.pop_front();
    }
    if (startNode == -1) {
        // can't find any useful links - just take any one from boundary then
        startNode = boundNodes.front().nodeNum;
    }


    list <int> toCheck; //nodes which need to be checked next
    list <int> filledBound; //container for nodes in order
    filledBound.push_back(startNode);
    flags[startNode] = 4;

    //cout << "Initial setting up." << endl;
    //2. put just one of the adjacent nodes into the toCheck list
    for (int linkCnt = 0; linkCnt < nNodesLinked[startNode]; linkCnt++) {
        int linkNode = linkMat[(startNode)*nMaxLinkedNodes+linkCnt];
        if (nodeBlocks[linkNode] == blockNum && flags[linkNode] == 2) {
            toCheck.push_back(linkNode);
            flags[linkNode] = 3;
            break;
        }
    }

    //cout << "Do flood fill." << endl;
    //3. flood fill as far as possible
    int doNode;
    int lastBlockLinked = -1;
    for (int testCnt = 0; testCnt < nNodes; testCnt++) {//infinite loop; make it smaller than nNodes just in case of errors/bugs
        if (toCheck.size() == 0) {
            //can't go any further...
            if (filledBound.size() == boundNodes.size()) {
                //all done!
                break;
            } else {
                //still have some left
                //cout << "Flood fill has reached end. Need to find new starting node " << filledBound.size() << endl;
                //need to find the most suitable next node(s) and set doNode
                if (bothBounds.size() > 0) {
                    doNode = -1;
                    for (list <int>:: iterator p = bothBounds.begin();
                         p != bothBounds.end(); p++) {

                        if (flags[*p] != 2) { //this node has already been taken so should be removed
                            p = bothBounds.erase(p);
                            p--;
                            continue;
                        }

                        // try to find one which is adjacent to the same block as last time
                        for (int linkCnt = 0; linkCnt < nNodesLinked[*p]; linkCnt++) {
                            int linkNode = linkMat[(*p)*nMaxLinkedNodes+linkCnt];
                            if (nodeBlocks[linkNode] == lastBlockLinked) {
                                doNode = *p;
                                flags[doNode] = 3;
                                p = bothBounds.erase(p);
                                p--;
                                break;
                            }
                        }
                        if (doNode != -1) {
                            //cout << "Found starting node (attached on bothBound list): " << doNode << endl;
                            break;
                        }
                    }

                    if (doNode == -1) { // haven't found a suitable one so just take first one off the list
                        doNode = bothBounds.front();
                        bothBounds.pop_front();
                        flags[doNode] = 3;
                        //cout << "Found starting node (first one on bothBound list): " << doNode << endl;
                    }


                } else {
                    doNode = -1;
                    //loop through boundary nodes finding an unused one that is linked to lastBlockLinked
                    for (list <node>:: iterator p = boundNodes.begin();
                         p != boundNodes.end(); p++) {
                        if (flags[p->nodeNum] == 2)
                            for (int linkCnt = 0; linkCnt < nNodesLinked[p->nodeNum]; linkCnt++) {
                                int linkNode = linkMat[(p->nodeNum)*nMaxLinkedNodes+linkCnt];
                                if (nodeBlocks[linkNode] == lastBlockLinked) {
                                    doNode = p->nodeNum;
                                    flags[p->nodeNum] = 3;
                                    //cout << "Found starting node (linked to same block on main list): " << doNode << endl;
                                    break;
                                }
                            }
                        if (doNode != -1) {
                            break;
                        }
                    }
                    if (doNode == -1) {
                        //loop through boundary taking first unused one
                        for (list <node>:: iterator p = boundNodes.begin();
                             p != boundNodes.end(); p++) {
                            if (flags[p->nodeNum] == 2) {
                                doNode = p->nodeNum;
                                //cout << "Found starting node (any on main list): " << doNode << endl;
                                flags[p->nodeNum] = 3;
                                break;
                            }
                        }
                        if (doNode == -1) {
                            stringstream errSs;

                            errSs << "Can't find free node! Help!" << endl;
                            int cnt = 0;
                            for (list<int>::iterator p = filledBound.begin();
                                 p != filledBound.end(); p++)
                                errSs << ++cnt << "," << *p << endl;

                            errSs << "Ordered size: " << filledBound.size() << ", original size: " << boundNodes.size() << endl;

                            err newError("Blocker - memory arrangement", 2411, errSs.str());
                            throw newError;
                        }
                    }
                }
                filledBound.push_back(doNode);
            }

        } else {
            //take it off the check list
            doNode = toCheck.front();
            toCheck.pop_front();

            //put it on the final list
            filledBound.push_back(doNode);
        }




        for (int linkCnt = 0; linkCnt < nNodesLinked[doNode]; linkCnt++) {
            int linkNode = linkMat[(doNode)*nMaxLinkedNodes+linkCnt];
            if (nodeBlocks[linkNode] == blockNum && flags[linkNode] == 2) {
                toCheck.push_back(linkNode);
                flags[linkNode] = 3;
            } else if (nodeBlocks[linkNode] != blockNum) {
                lastBlockLinked = blockNum; //store this for later
            }
        }
        flags[doNode] = 4;
    }


    //5. move loop around to get good starting point

    // get a list of surrounding blocks
    list<surrBlockCnt> blockNums;
    for (list<int> :: iterator p = filledBound.begin();
         p != filledBound.end();
         p++) {
        for (int linkCnt = 0; linkCnt < nNodesLinked[*p]; linkCnt++) {
            int linkNode = linkMat[(*p)*nMaxLinkedNodes+linkCnt];
            int linkedBlock = nodeBlocks[linkNode];
            if (linkedBlock == blockNum)
                continue;
            bool alreadyExists = false;
            if (blockNums.size() > 0)
                for (list<surrBlockCnt> :: iterator q = blockNums.begin();
                     q != blockNums.end();
                     q++) {
                    if (q->blockNum == linkedBlock) {
                        alreadyExists = true;
                    }
                }
            if (!alreadyExists) {
                surrBlockCnt newSurr;
                newSurr.blockNum = linkedBlock;
                //newSurr.cnt++;
                blockNums.push_back(newSurr);
            }
        }
    }
    // total up how many times each block is linked to
    for (list<int> :: iterator p = filledBound.begin();
         p != filledBound.end();
         p++) {
        for (list<surrBlockCnt> :: iterator q = blockNums.begin();
             q != blockNums.end();
             q++) {
            for (int linkCnt = 0; linkCnt < nNodesLinked[*p]; linkCnt++) {
                int linkNode = linkMat[(*p)*nMaxLinkedNodes+linkCnt];
                int linkedBlock = nodeBlocks[linkNode];
                if (q->blockNum == linkedBlock) {
                    (q->cnt)++;
                    break;
                }
            }

        }
    }
    //find the maximum
    int maxSize = 0;
    int maxBlock = -1;
    for (list<surrBlockCnt> :: iterator q = blockNums.begin();
         q != blockNums.end();
         q++) {
        if (q->cnt > maxSize) {
            maxSize = q->cnt;
            maxBlock = q->blockNum;
        }
    }


    //6. put in order so block with most links comes first

    //need to find start of maxBlock.
    //check end of list to see whether that is 'part' of maxBlock
    list<int> :: iterator p = filledBound.end();
    p--;
    bool foundEnd = false;
    for (int linkCnt = 0; linkCnt < nNodesLinked[*p]; linkCnt++) {
        int linkNode = linkMat[(*p)*nMaxLinkedNodes+linkCnt];
        int linkedBlock = nodeBlocks[linkNode-1];
        if (linkedBlock == maxBlock) {
            //we have found it!
            foundEnd = true;
            break;
        }
    }

    if (foundEnd == true) { // if so we join accordingly
        p = filledBound.end();
        p--;
        for (; p != filledBound.begin(); p--) {//loop through to the beginning of the block
            bool found = false;
            for (int linkCnt = 0; linkCnt < nNodesLinked[*p]; linkCnt++) {
                int linkNode = linkMat[(*p)*nMaxLinkedNodes+linkCnt];
                int linkedBlock = nodeBlocks[linkNode];
                if (linkedBlock == maxBlock) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                //this is the first node where maxBlock hasn't been linked to
                p++; //increase by 1 so that p links to the first one
                break;
            }
        }
    } else {
        p = filledBound.begin();
        for (; p != filledBound.end(); p++) {
            bool found = false;
            for (int linkCnt = 0; linkCnt < nNodesLinked[*p]; linkCnt++) {
                int linkNode = linkMat[(*p)*nMaxLinkedNodes+linkCnt];
                int linkedBlock = nodeBlocks[linkNode];
                if (linkedBlock == maxBlock) {
                    //we have found it!
                    found = true;
                    break;
                }
            }
            if (found) break;
        }
    }


    //so now, p points to the first node attached to maxBlock.
    list<int> boundOrdered;
    boundOrdered.splice(boundOrdered.begin(), filledBound, p, filledBound.end());
    boundOrdered.splice(boundOrdered.end(), filledBound);


    //7. Put into block data structure
    for (int cnt = 0; cnt < XBLOCKSIZE*YBLOCKSIZE; cnt++)
        nodeRef[cnt] = -1;


    //do linked data first
    p = boundOrdered.begin();
    for (unsigned int cnt = 0; cnt < boundOrdered.size(); cnt++) {
        if (*p > nNodes) {
            stringstream errSs;
            errSs << "Error - referenced node too large" << endl;
            err newError("Blocker - memory arrangement", 2412, errSs.str());
            throw newError;
        }
        nodeRef[cnt] = *p;
        p++;
        if (p == boundOrdered.end()) break;
    }

    //do internal only data
    if (intNodes.size() > 0) {
        list <node>:: iterator p = intNodes.begin();
        for (int cnt = boundOrdered.size(); cnt < XBLOCKSIZE*YBLOCKSIZE; cnt++) {
            nodeRef[cnt] = p->nodeNum;
            p++;
            if (p == intNodes.end()) break;
        }
    }

    //blocksList.push_back(newBlockData);
    blockFlag = 2;


    return 0;
}

//-------------------------------------------------------------------------------------------------------------------

int block::updateColsRows(int * cols, int * rows) {
    if (blockFlag < 2) {
        stringstream errSs;
        errSs << "Error - no node refs defined " << endl;
        err newError("Blocker - updateColsRows", 2412, errSs.str());
        throw newError;
    }
    for (int colCnt = 0; colCnt < XBLOCKSIZE; colCnt++)
        for (int rowCnt = 0; rowCnt < YBLOCKSIZE; rowCnt++) {
            int node = nodeRef[colCnt + rowCnt*XBLOCKSIZE];
            //blockMap[colCnt + rowCnt*32 + blockCnt*32*16] = node;
            if (node >= 0) {
                //cout << node << ":" << colCnt << "," << rowCnt << endl;
                if (node >= nNodes) {
                    stringstream errSs;
                    errSs << "Error in node number: " << node << endl;
                    errSs << blockNum << "," << colCnt << "," << rowCnt << endl;
                    err newError("Blocker - updateColsRows", 2413, errSs.str());
                    throw newError;
                }
                cols[node] = colCnt;
                rows[node] = rowCnt;
            }
        }
    return 0;
}

//-------------------------------------------------------------------------------------------------------------------

int block::allocateLoads(int cnt, int * loadMaps, int nLinkedBlocks, int * connArray) {
    if (cnt >= NBLOCKLINKS) {
        return -1;
    }
    int loadLength = loadLengths[cnt];
    //number of positions within a row the load can be allocated to
    int nRowPos = NACCESSCOLS-loadLength+1;

    int maxPriRow = -1;
    int maxPriCol = -1;
    int maxPriVal = -1;

    int checkMask = (1 << loadLength)-1;

    for (int rowCnt = 0; rowCnt < nLinkedBlocks*NACCESSROWS; rowCnt++) {
        for (int colCnt = 0; colCnt < nRowPos; colCnt++) {
            //check the fitting
            int binComb = ((connArray[rowCnt] >> (nRowPos-1-colCnt)) & checkMask);
            int totBits = 0;
            for (int bitCnt = 0; bitCnt < loadLength; bitCnt++) {
                if (binComb & 1) totBits++;
                binComb >>= 1;
            }
            if (totBits > maxPriVal) {
                maxPriVal = totBits;
                maxPriCol = colCnt;
                maxPriRow = rowCnt;
            }
        }
    }
    //allocate it to the highest priority one
    int rowBefore = connArray[maxPriRow];
    connArray[maxPriRow] -= connArray[maxPriRow] & (checkMask << (nRowPos-1-maxPriCol));
    //check whether they are all zeros
    bool allzeros = true;
    for (int rowCnt = 0; rowCnt < nLinkedBlocks*NACCESSROWS; rowCnt++) {
        if (connArray[rowCnt] != 0) {
            allzeros = false;
            break;
        }
    }
    if (allzeros == true) {
        loadMaps[cnt] = maxPriRow*NACCESSCOLS + (nRowPos-1-maxPriCol);
        return 0;
    }
    int res = allocateLoads(cnt+1, loadMaps, nLinkedBlocks, connArray);
    if (res == -1) {
        connArray[maxPriRow] = rowBefore;
        return -1; //could possibly try allocating another way rather than just giving up
    }
    if (res == 0) {
        loadMaps[cnt] = maxPriRow*NACCESSCOLS + (nRowPos-1-maxPriCol);
        return 0;
    }
    return 2;
}

//-------------------------------------------------------------------------------------------------------------------


int block::linkArrange(int * cols, int * rows, vector<block> &blocks) {
    //return -10 means we can improve things by splitting
    //return -2 means there are problems with the linked blocks (i.e. not ready yet)

    if (blockFlag >= 3) {
        // don't do anything because block is already finalised
        return 0;
    }

    if (blockFlag < 2) {
        stringstream errSs;
        errSs << "No node refs defined " << blockNum << endl;
        err newError("Blocker - linkArrange", 2414, errSs.str());
        throw newError;
    }

    //get the adjacent blocks
    //need to do this to check that all surrounding blocks are finalised
    list <int> attachedBlocks;
    getLinks(attachedBlocks);
    for (list<int>::iterator p = attachedBlocks.begin();
         p != attachedBlocks.end(); p++) {
        if (blocks[*p].blockFlag < 2) {
            return -2;
        }
    }

    for (int cnt = 0; cnt < NBLOCKLINKS; cnt++) {
        blockPointers[cnt] = 0;
    }


    list<int> linkedNodes;
    vector<blockLink> linkedBlocks;
    //NB:
//    struct blockLink {
//        int blockNum;
//        list <int> nodes;
//    };
    // this stores the nodes for each linked block

    for (list<node>:: iterator n = boundNodes.begin(); n!=boundNodes.end(); n++) {
        int cnt = n->nodeNum;
            for (int linkCnt = 0; linkCnt < nNodesLinked[cnt]; linkCnt++) {
                int linkNode = linkMat[cnt*nMaxLinkedNodes+linkCnt];
                if (nodeBlocks[linkNode] == blockNum)
                    continue;

                //here it must come from another block

                //search through and see if this node is already in the list
                bool alreadyListed = false;
                if (linkedNodes.size() > 0)
                    for (list<int>::iterator p = linkedNodes.begin();
                         p!=linkedNodes.end(); p++) {
                        if (linkNode == *p) {
                            alreadyListed = true;
                            break;
                        }
                    }

                //if it hasn't already been put into the list, add it to the end
                if (!alreadyListed) {
                    linkedNodes.push_back(linkNode);

                    //add the block to the block list if not already listed
                    int linkBlock = nodeBlocks[linkNode];
                    //search through and see if it's already in the list
                    bool alreadyListedBlock = false;
                    if (linkedBlocks.size() > 0)
                        for (vector<blockLink>::iterator q = linkedBlocks.begin();
                             q!=linkedBlocks.end(); q++) {
                            if (linkBlock == q->blockNum) {
                                alreadyListedBlock = true;
                                q->nodes.push_back(linkNode);
                                break;
                            }
                        }
                    if (!alreadyListedBlock) {
                        blockLink addBlock;
                        addBlock.blockNum = linkBlock;
                        addBlock.nodes.push_back(linkNode);
                        linkedBlocks.push_back(addBlock);
                    }
                }
            }
    }
    if (linkedNodes.size() == 0) {
        stringstream errSs;
        errSs << "No linked nodes found for " << blockNum << endl;
        err newError("Blocker - linkArrange", 2415, errSs.str());
        throw newError;
    }





    int nLinkedBlocks = linkedBlocks.size();
    int * connArray = new int[NACCESSROWS*nLinkedBlocks];

    for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++)
        connArray[cnt] = 0;

    //build up binary values for each row
    for (int cnt = 0; cnt < nLinkedBlocks; cnt++) {
        for (list<int>::iterator p = linkedBlocks[cnt].nodes.begin();
             p != linkedBlocks[cnt].nodes.end(); p++) {
            if (rows[*p] > NACCESSROWS) {//trying to link to another block beyond the specified number of rows
                return -10;
            }
            connArray[rows[*p] + cnt*NACCESSROWS] |= 1<<(cols[*p]>>3);
        }
    }




    int loadPos = 0;

    int loadMaps[NBLOCKLINKS];
    for (int cnt = 0; cnt < NBLOCKLINKS; cnt++)
        loadMaps[cnt] = 0;


    while (loadPos < NBLOCKLINKS) {
        if (loadLengths[loadPos] == 4) {
            //check completely filled rows
            for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++)
                if (connArray[cnt] == 15) {
                    loadMaps[loadPos] = cnt*4+0;
                    loadPos++;
                    connArray[cnt] = 0;
                    if (loadLengths[loadPos] != 4)
                        break;
                }
        }
        if (loadLengths[loadPos] == 4) {
            //check 3/4 filled rows
            for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++)
                if (connArray[cnt] == 14 || connArray[cnt] == 13 || connArray[cnt] == 11 || connArray[cnt] == 7) {
                    loadMaps[loadPos] = cnt*4+0;
                    loadPos++;
                    connArray[cnt] = 0;
                    if (loadLengths[loadPos] != 4)
                        break;
                }
        }
        if (loadLengths[loadPos] == 4) {
            //check 2/4, split
            for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++)
                if (connArray[cnt] == 9 || connArray[cnt] == 10 || connArray[cnt] == 5) {
                    loadMaps[loadPos] = cnt*4+0;
                    loadPos++;
                    connArray[cnt] = 0;
                    if (loadLengths[loadPos] != 4)
                        break;
                }
        }
        if (loadLengths[loadPos] == 4) {
            //check 2/4, joined
            for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++)
                if (connArray[cnt] == 3 || connArray[cnt] == 6 || connArray[cnt] == 12) {
                    loadMaps[loadPos] = cnt*4+0;
                    loadPos++;
                    connArray[cnt] = 0;
                    if (loadLengths[loadPos] != 4)
                        break;
                }
        }
        if (loadLengths[loadPos] == 4) {
            //check 1/4
            for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++)
                if (connArray[cnt] == 1 || connArray[cnt] == 2 || connArray[cnt] == 4 || connArray[cnt] == 8) {
                    loadMaps[loadPos] = cnt*4+0;
                    loadPos++;
                    connArray[cnt] = 0;
                    if (loadLengths[loadPos] != 4)
                        break;
                }
        }

        if (loadLengths[loadPos] == 2) {
            //check 2/4, joined
            for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++) {
                if ((connArray[cnt]&3) == 3) {
                    loadMaps[loadPos] = cnt*4+0;
                    loadPos++;
                    connArray[cnt] &= 12;
                    if (loadLengths[loadPos] != 2)
                        break;
                }
                if ((connArray[cnt]&6) == 6) {
                    loadMaps[loadPos] = cnt*4+1;
                    loadPos++;
                    connArray[cnt] &= 9;
                    if (loadLengths[loadPos] != 2)
                        break;
                }

                if ((connArray[cnt]&12) == 12) {
                    loadMaps[loadPos] = cnt*4+2;
                    loadPos++;
                    connArray[cnt] &= 3;
                    if (loadLengths[loadPos] != 2)
                        break;
                }

            }
        }
        if (loadLengths[loadPos] == 2) {
            //check 1/4
            for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++) {
                if ((connArray[cnt]&1) == 1) {
                    loadMaps[loadPos] = cnt*4+0;
                    loadPos++;
                    connArray[cnt] &= 12;
                    if (loadLengths[loadPos] != 2)
                        break;
                }
                if ((connArray[cnt]&2) == 2) {
                    loadMaps[loadPos] = cnt*4+0;
                    loadPos++;
                    connArray[cnt] &= 12;
                    if (loadLengths[loadPos] != 2)
                        break;
                }

                if ((connArray[cnt]&4) == 4) {
                    loadMaps[loadPos] = cnt*4+2;
                    loadPos++;
                    connArray[cnt] &= 3;
                    if (loadLengths[loadPos] != 2)
                        break;
                }
                if ((connArray[cnt]&8) == 8) {
                    loadMaps[loadPos] = cnt*4+2;
                    loadPos++;
                    connArray[cnt] &= 3;
                    if (loadLengths[loadPos] != 2)
                        break;
                }
            }
        }
        if (loadLengths[loadPos] == 1) {
            //check 1/4
            for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++) {
                if (connArray[cnt] != 0) {
                    for (int cntBit = 0; cntBit < 4; cntBit++) {
                        if (connArray[cnt] & (1 << cntBit)) {
                            loadMaps[loadPos] = cnt*4+cntBit;
                            loadPos++;
                            connArray[cnt] &= 15-(1<<cntBit);
                            if (loadPos >= 16)
                                break;
                        }
                    }
                    if (loadPos >= NBLOCKLINKS)
                        break;
                }
            }
        }


        bool done = true;
        for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++)
            if (connArray[cnt] != 0) {
                done = false;
                break;
            }
        if (done == true)
            break;

        if (loadPos >= 16) {
            cout << "Too many loads required from block " << blockNum << ". Rearranging." << endl;
            return -10;
        }
    }

    //check to see if we've used up all the memory
    bool done = true;
    for (int cnt = 0; cnt < nLinkedBlocks*NACCESSROWS; cnt++)
        if (connArray[cnt] != 0) {
            done = false;
            break;
        }
    if (done == false) {
        cout << "Problem. Not enough room." << endl;
        return -10;
        //return;
    }


//    int res = allocateLoads(0,loadMaps,nLinkedBlocks,connArray);
//    if (res == -1) {
//        cout << "Too many loads required from block " << blockNum << ". Rearranging." << endl;
//        return -10;
//    }





    //now need to put it into blocks

    for (int cnt = 0; cnt < NBLOCKLINKS; cnt++) {
        int totRow = loadMaps[cnt] / NACCESSCOLS;
        int locBlock = totRow / NACCESSROWS;
        int row = totRow % NACCESSROWS;
        int globBlockNum = linkedBlocks[locBlock].blockNum;
        int col = loadMaps[cnt] % NACCESSCOLS;

        blockPointers[cnt] = (globBlockNum << 8) + (row << 4) + col;
    }

    blockFlag = 3;



    delete [] connArray;
    return 0;
}

void block::removeOrderFlags() {
    blockFlag = 1;

    for (list<node>::iterator p = intNodes.begin();
         p!= intNodes.end(); p++) {
        flags[p->nodeNum] = 1;
    }

    for (list<node>::iterator p = boundNodes.begin();
         p!= boundNodes.end(); p++) {
        flags[p->nodeNum] = 2;
    }

    for (int cnt = 0; cnt < XBLOCKSIZE*YBLOCKSIZE; cnt++)
        nodeRef[cnt] = -1;

    for (int cnt = 0; cnt < NBLOCKLINKS; cnt++)
        blockPointers[cnt] = 0;

}

void block::removeLinkFlags() {
    if (blockFlag <= 2)
        return;

    blockFlag = 2;

//    for (list<node>::iterator p = intNodes.begin();
//         p!= intNodes.end(); p++) {
//        flags[p->nodeNum] = 1;
//    }

//    for (list<node>::iterator p = boundNodes.begin();
//         p!= boundNodes.end(); p++) {
//        flags[p->nodeNum] = 2;
//    }

    for (int cnt = 0; cnt < NBLOCKLINKS; cnt++)
        blockPointers[cnt] = 0;
}
