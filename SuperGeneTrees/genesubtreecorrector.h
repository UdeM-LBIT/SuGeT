#ifndef GENESUBTREECORRECTOR_H
#define GENESUBTREECORRECTOR_H

#include "trees/node.h"
#include "trees/genespeciestreeutil.h"
#include "supergenetreemaker.h"


#include <unordered_map>
#include <vector>
#include <algorithm>

using namespace std;

class GeneSubtreeCorrector
{
public:
    GeneSubtreeCorrector();

    Node* GetSubtreeTripletRespectingHistory(Node* geneTree, Node* speciesTree,
                                             unordered_map<Node*, Node*> lca_mapping,
                                             vector<Node *> geneSubtreeRoots, bool mustRetainDupSpec);
};

#endif // GENESUBTREECORRECTOR_H