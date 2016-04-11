
#include "trees/node.h"
#include "trees/newicklex.h"
#include "trees/genespeciestreeutil.h"
#include "supergenetreemaker.h"
#include "genesubtreecorrector.h"

#include <iostream>
#include <string>
#include <map>

using namespace std;


void DoSuperGeneTree(string gcontent, string scontent, bool preserveDupSpec);
void DoSubtreeCorrection(string gcontent, string scontent, bool preserveDupSpec);

int main(int argc, char *argv[])
{
    /*
    here I store my previous test command line args
    -p -gn "((A1__A, C1__C),(B1__B, (C2__C,D2__D)));((B2__B,D3__D),(B1__B,C3__C));" -sn "((A,B),(C,D))"
    -p -m sub -gn "(((A1__A, B1__B)m, (A2__A, B2__B)m), ((C3__C, D3__D), ((A4__A,B4__B),(C4__C,D4__D)))m);" -sn "((A,B),(C,D))"

    -m sub -sn "(((A,B),C), (D,E));" -gn "(((A1__A, B1__B)m,C1__Cm), ((D1__D,E1__E),(((A2__A,B2__B),C2__C), (D2__D,E2__E)))m);"

    -m sub -sn "(((A,B),C), (D,E));" -gn "(((((A1__A, D1__D),B1__B)m,((B2__B, C2__C),E2__E)m),(((D3__D, E3__E),(D4__D, E4__E))m,((A5__A, B5__B), (C5__C, C6__C))m)),((A6__A, D6__D),(B7__B, (B8__B, C8__C)))m);"

    */

    //gene tree label format is GENENAME__SPECIESNAME

    //below = supertree mode
    int verbose = 0;
    bool preserveDupSpec = false;
    map<string, string> args;

    string prevArg = "";
    for (int i = 0; i < argc; i++)
    {
        if (string(argv[i]) == "-v")
        {
            verbose = 1;
            prevArg = "";
        }
        else if (string(argv[i]) == "-v2")
        {
            verbose = 2;
            prevArg = "";
        }
        else if (string(argv[i]) == "-p")   //PRESERVE
        {
            preserveDupSpec = true;
            prevArg = "";
        }
        else
        {
            if (prevArg != "" && prevArg[0] == '-')
            {
                args[Util::ReplaceAll(prevArg, "-", "")] = string(argv[i]);
            }

            prevArg = string(argv[i]);
        }
    }


    string mode = "sgt";

    if (args.find("m") != args.end())
    {
        mode = args["m"];

        if (mode != "sgt" && mode != "sub")
        {
            cout<<"Mode must be one of 'sgt' or 'sub'"<<endl;
            return 0;
        }
    }


    string scontent = "";
    string gcontent = "";

    if (args.find("s") != args.end())
    {
        scontent = Util::GetFileContent(args["s"]);
    }
    else if (args.find("sn") != args.end())
    {
        scontent   = args["sn"];
    }



    if (args.find("g") != args.end())
    {
        gcontent = Util::GetFileContent(args["g"]);
    }
    else if (args.find("gn") != args.end())
    {
        gcontent = args["gn"];
    }


    if (mode == "sgt")
    {
        DoSuperGeneTree(gcontent, scontent, preserveDupSpec);
    }
    else if (mode == "sub")
    {
        DoSubtreeCorrection(gcontent, scontent, preserveDupSpec);
    }


}



void DoSubtreeCorrection(string gcontent, string scontent, bool preserveDupSpec)
{
    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);

    Node* geneTree = NewickLex::ParseNewickString(gcontent, false);

    //find subtree roots.  Their label is 'm'
    vector<Node*> markedGeneTreeNodes;
    TreeIterator* it = geneTree->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        string lbl = n->GetLabel();
        if (lbl.length() > 0 && Util::Trim(lbl).substr(lbl.length() - 1, 1) == "m")
        {
            markedGeneTreeNodes.push_back(n);

            //remove extra m, so we don't mess gene species mapping
            if (n->IsLeaf())
                n->SetLabel(lbl.substr(0, lbl.length() - 1));
        }
    }
    geneTree->CloseIterator(it);

    unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);


    GeneSubtreeCorrector gsc;
    Node* newgenetree = gsc.GetSubtreeTripletRespectingHistory(geneTree, speciesTree, lcamap, markedGeneTreeNodes, preserveDupSpec);



    if (!newgenetree)
    {
        cout<<"No tree was returned.  This should never happen.  Never.  The fact that you are facing the physically impossible makes us question the laws of the universe which were believed true."<<endl;
    }
    else
    {
        cout<<NewickLex::ToNewickString(newgenetree)<<endl;
        delete newgenetree;
    }

    delete geneTree;
    delete speciesTree;
}


void DoSuperGeneTree(string gcontent, string scontent, bool preserveDupSpec)
{
    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);

    vector<string> strGeneTrees = Util::Split(gcontent, ";", false);
    vector<Node*> geneTrees;
    vector< unordered_map<Node*, Node*> > lca_mappings;

    for (int i = 0; i < strGeneTrees.size(); i++)
    {
        string tmpnewick = Util::ReplaceAll(strGeneTrees[i], "\n", "");
        Node* geneTree = NewickLex::ParseNewickString(tmpnewick);

        unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);


        geneTrees.push_back(geneTree);
        lca_mappings.push_back(lcamap);


        TreeIterator* it = geneTree->GetPostOrderIterator();
        while (Node* n = it->next())
        {
            Node* smap = lcamap[n];

            if (!smap)
            {
                cout<<n->GetLabel()<<" is not mapped!";
                return;
            }
        }
        geneTree->CloseIterator(it);
    }


    SuperGeneTreeMaker sgtMaker;
    pair<Node*, int> res = sgtMaker.GetSuperGeneTreeMinDL(geneTrees, lca_mappings, speciesTree, preserveDupSpec, true);

    Node* superTree = res.first;
    int cost = res.second;


    if (!superTree)
    {
        cout<<"It seems that no solution exists."<<endl;
    }
    else
    {
        cout<<"DLCOST="<<cost<<endl;
        cout<<NewickLex::ToNewickString(superTree)<<endl;

        delete superTree;
    }

    for (int i = 0; i < geneTrees.size(); i++)
        delete geneTrees[i];
    delete speciesTree;
}
