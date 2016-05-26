
#include "trees/node.h"
#include "trees/newicklex.h"
#include "trees/genespeciestreeutil.h"
#include "supergenetreemaker.h"
#include "genesubtreecorrector.h"
#include "trees/polysolver.h"

#include <iostream>
#include <string>
#include <map>
#include <time.h>
#include "div/tinydir.h"



using namespace std;


string DoSuperGeneTree(string gcontent, string scontent, bool preserveDupSpec, string outputmode = "");
string DoSubtreeCorrection(string gcontent, string scontent, bool preserveDupSpec, string markedNodesMode = "", string outputmode = "tree");


string DoSGTOnHighSpecs(string gcontent, string scontent, bool preserveDupSpec);
void DoPolytomyCorrection(string gcontent, string scontent);


int main(int argc, char *argv[])
{
    /*
    here I store my previous test command line args
    -p -gn "((A1__A, C1__C),(B1__B, (C2__C,D2__D)));((B2__B,D3__D),(B1__B,C3__C));" -sn "((A,B),(C,D))"
    -p -m sub -gn "(((A1__A, B1__B)m, (A2__A, B2__B)m), ((C3__C, D3__D), ((A4__A,B4__B),(C4__C,D4__D)))m);" -sn "((A,B),(C,D))"

    -m sub -sn "(((A,B),C), (D,E));" -gn "(((A1__A, B1__B)m,C1__Cm), ((D1__D,E1__E),(((A2__A,B2__B),C2__C), (D2__D,E2__E)))m);"

    -m sub -sn "(((A,B),C), (D,E));" -gn "(((((A1__A, D1__D),B1__B)m,((B2__B, C2__C),E2__E)m),(((D3__D, E3__E),(D4__D, E4__E))m,((A5__A, B5__B), (C5__C, C6__C))m)),((A6__A, D6__D),(B7__B, (B8__B, C8__C)))m);"

    BATCH
    -m sub -s "C:/Users/Manuel/Desktop/projects/code/data/species_tree.ensembl.topology.nw" -b "C:/Users/Manuel/Desktop/projects/code/data/gene_trees_may2016/"

    ENSGT00460000041710.txt
    -m sub -p -s "C:/Users/Manuel/Desktop/projects/code/data/species_tree.ensembl.topology.nw" -g "C:/Users/Manuel/Desktop/projects/code/data/gene_trees_may2016/ENSGT00460000041710.txt"


    FUCK THIS BELOW
    -m sub -s "C:/Users/Manuel/Desktop/projects/code/data/species_tree.ensembl.topology.nw" -b "C:/Users/Manuel/Desktop/projects/code/data/test_gene_trees/"

    -m sub -s "C:/Users/Manuel/Desktop/projects/code/data/species_tree.ensembl.topology.nw" -b "C:/Users/Manuel/Desktop/projects/code/data/test_gene_trees/"
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

    string batchDir = "";
    if (args.find("b") != args.end())  //BATCH DIR
    {
        batchDir = args["b"];
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


    //////////////////////////////////////////////////////////
    // BATCH MODE
    //////////////////////////////////////////////////////////
    if (batchDir != "")
    {
        batchDir = Util::ReplaceAll(batchDir, "\\", "/");
        batchDir += "/";

        tinydir_dir dir;
        tinydir_open(&dir, batchDir.c_str());

        int cpt = 0;

        //header line
        cout<<"FILE,"
            <<"NB_leaves,"
            <<"NB_marked_nodes,"
            <<"DUPS_before,"
            <<"LOSSES_before,"
            <<"RECONC_before,"
            <<"DUPS_after,"
            <<"LOSSES_after,"
            <<"RECONC_after,"
            <<"NB_marked_after,"
            <<"TIME(ms),"
            <<"NB_leaves_lbl,"
            <<"NB_marked_nodes_lbl,"
            <<"DUPS_before_lbl,"
            <<"LOSSES_before_lbl,"
            <<"RECONC_before_lbl,"
            <<"DUPS_after_lbl,"
            <<"LOSSES_after_lbl,"
            <<"RECONC_after_lbl,"
            <<"NB_marked_after_lbl,"
            <<"TIME_lbl,"
            <<"Poly_DUPS,"
            <<"Poly_LOSSES,"
            <<"Poly_REC,"
            <<"Poly_TIME,"
            //<<"SGT_DUPS,"
            //<<"SGT_LOSSES,"
            //<<"SGT_REC,"
            //<<"SGT_TIME,"
            <<"SGT_DUPS_lbl,"
            <<"SGT_LOSSES_lbl,"
            <<"SGT_REC_lbl,"
            <<"SGT_TIME_lbl"
            <<endl;

        while (dir.has_next)
        {
            tinydir_file file;
            tinydir_readfile(&dir, &file);


            if (file.is_dir)
            {
                ;
            }
            else
            {

                string filename = string(file.name);    //for some reason, calling file.name later causes crashes
                string fullname = batchDir + filename;


                //cout<<"Parsing "<<fullname<<endl;

                gcontent = Util::GetFileContent(fullname);

                string pruned_scontent = GeneSpeciesTreeUtil::Instance()->GetPrunedSpeciesTreeNewick(gcontent, scontent);

                if (mode == "sgt")
                {
                    DoSuperGeneTree(gcontent, pruned_scontent, preserveDupSpec);
                }
                else if (mode == "sub")
                {
                    cout<<filename<<",";

                    ///////////////////////////////////////////////////////////
                    /// ROUND 1 : MTRH without labels
                    ///////////////////////////////////////////////////////////
                    clock_t time = clock();
                    DoSubtreeCorrection(gcontent, pruned_scontent, false, "highspecs", "stats");

                    time = clock() - time;
                    int ms = (double)time / CLOCKS_PER_SEC * 1000;
                    cout<<","<<ms<<",";

                    ///////////////////////////////////////////////////////////
                    /// ROUND 2 : MTRH with labels
                    ///////////////////////////////////////////////////////////
                    time = clock();
                    DoSubtreeCorrection(gcontent, pruned_scontent, true, "highspecs", "stats");  //repeating useless preprocessing here but do I care?

                    time = clock() - time;
                    ms = (double)time / CLOCKS_PER_SEC * 1000;
                    cout<<","<<ms<<",";


                    ///////////////////////////////////////////////////////////
                    /// ROUND 3 : PolytomySolver, for Nadia
                    ///////////////////////////////////////////////////////////
                    time = clock();

                    DoPolytomyCorrection(gcontent, pruned_scontent);

                    time = clock() - time;
                    ms = (double)time / CLOCKS_PER_SEC * 1000;
                    cout<<","<<ms<<",";

                    ///////////////////////////////////////////////////////////
                    /// ROUND 4 : SuGeT  without labels
                    ///////////////////////////////////////////////////////////
                    /// NOPE, TOO LONG
                    time = clock();

                    DoSGTOnHighSpecs(gcontent, pruned_scontent, false);

                    time = clock() - time;
                    ms = (double)time / CLOCKS_PER_SEC * 1000;
                    cout<<","<<ms<<",";
                    ///////////////////////////////////////////////////////////
                    /// ROUND 5 : SuGeT  with labels
                    ///////////////////////////////////////////////////////////
                    time = clock();

                    DoSGTOnHighSpecs(gcontent, pruned_scontent, true);

                    time = clock() - time;
                    ms = (double)time / CLOCKS_PER_SEC * 1000;
                    cout<<","<<ms;

                    cout<<endl;
                }
            }

            cpt++;


            tinydir_next(&dir);
        }

        tinydir_close(&dir);
    }
    //////////////////////////////////////////////////////////
    // SINGLE GENE TREE MODE
    //////////////////////////////////////////////////////////
    else
    {
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


}



string DoSubtreeCorrection(string gcontent, string scontent, bool preserveDupSpec, string markedNodesMode, string outputmode)
{

    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);
    Node* geneTree = NewickLex::ParseNewickString(gcontent, false);


    if (outputmode == "stats")
        cout<<geneTree->GetLeafSet().size()<<",";

    unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);

    GeneSpeciesTreeUtil::Instance()->GetDLScore(geneTree, speciesTree, lcamap);

    int DUPS_before = GeneSpeciesTreeUtil::Instance()->LASTNBDUPS;
    int LOSSES_before = GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES;





    //find subtree roots.  Their label is 'm', unless markedNodesMode = 'highspecs'
    vector<Node*> markedGeneTreeNodes;

    if (markedNodesMode == "")
    {
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
    }
    else if (markedNodesMode == "highspecs")
    {
        markedGeneTreeNodes = GeneSpeciesTreeUtil::Instance()->GetGeneTreeHighestSpeciations(geneTree, speciesTree, lcamap);
    }


    if (outputmode == "stats")
    {
        cout<<markedGeneTreeNodes.size()<<",";
    }


    GeneSubtreeCorrector gsc;
    Node* newgenetree = gsc.GetSubtreeTripletRespectingHistory(geneTree, speciesTree, lcamap, markedGeneTreeNodes, preserveDupSpec);

    string result = "";

    if (!newgenetree)
    {
        cout<<"No tree was returned.  This should never happen.  Never.  The fact that you are facing the physically impossible makes us question the laws of the universe which were believed true."<<endl;
    }
    else
    {
        result = NewickLex::ToNewickString(newgenetree);

        //one last reconciliation cost
        unordered_map<Node*, Node*> lcamap_new = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(newgenetree, speciesTree, "__", 1);
        GeneSpeciesTreeUtil::Instance()->GetDLScore(newgenetree, speciesTree, lcamap_new);
        int DUPS_aftercorrection = GeneSpeciesTreeUtil::Instance()->LASTNBDUPS;
        int LOSSES_aftercorrection = GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES;



        if ( outputmode == "stats" )
        {
                vector<Node*> marked_new = GeneSpeciesTreeUtil::Instance()->GetGeneTreeHighestSpeciations(newgenetree, speciesTree, lcamap_new);

                cout<<DUPS_before<<","
                    <<LOSSES_before<<","
                    <<(DUPS_before + LOSSES_before)<<","
                    <<DUPS_aftercorrection<<","
                    <<LOSSES_aftercorrection<<","
                    <<(DUPS_aftercorrection + LOSSES_aftercorrection)<<","
                    <<marked_new.size();

        }
        else
        {

            cout<<result<<endl;
        }

        //cout<<NewickLex::ToNewickString(newgenetree)<<endl;
        delete newgenetree;
    }

    delete geneTree;
    delete speciesTree;

    return result;

}


string DoSuperGeneTree(string gcontent, string scontent, bool preserveDupSpec, string outputmode)
{
    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);

    vector<string> strGeneTrees = Util::Split(gcontent, ";", false);
    vector<Node*> geneTrees;
    vector< unordered_map<Node*, Node*> > lca_mappings;

    for (int i = 0; i < strGeneTrees.size(); i++)
    {
        string tmpnewick = Util::ReplaceAll(strGeneTrees[i], "\n", "") + ";";
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
                return "";
            }
        }
        geneTree->CloseIterator(it);
    }


    SuperGeneTreeMaker sgtMaker;
    pair<Node*, int> res = sgtMaker.GetSuperGeneTreeMinDL(geneTrees, lca_mappings, speciesTree, preserveDupSpec, true);

    Node* superTree = res.first;
    int cost = res.second;


    string result = "";
    if (!superTree)
    {
        cout<<"It seems that no solution exists."<<endl;
    }
    else
    {
        result = NewickLex::ToNewickString(superTree);
        if (outputmode == "")
        {
            cout<<"DLCOST="<<cost<<endl;
            cout<<result<<endl;
        }
        else if (outputmode == "stats")
        {
            //stats: dup, loss, rec
            unordered_map<Node*, Node*> lcamap_new = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(superTree, speciesTree, "__", 1);
            int dl = GeneSpeciesTreeUtil::Instance()->GetDLScore(superTree, speciesTree, lcamap_new);

            cout<<GeneSpeciesTreeUtil::Instance()->LASTNBDUPS<<","
                <<GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES<<","
                <<dl;
        }

        delete superTree;
    }

    for (int i = 0; i < geneTrees.size(); i++)
        delete geneTrees[i];
    delete speciesTree;

    return result;
}





void DoPolytomyCorrection(string gcontent, string scontent)
{
    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);
    Node* geneTree = NewickLex::ParseNewickString(gcontent, false);


    unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);


    vector<Node*> highspecs = GeneSpeciesTreeUtil::Instance()->GetGeneTreeHighestSpeciations(geneTree, speciesTree, lcamap);

    if (highspecs.size() <= 2)
    {
        //nothing to do...
        int dl = GeneSpeciesTreeUtil::Instance()->GetDLScore(geneTree, speciesTree, lcamap);

        cout<<GeneSpeciesTreeUtil::Instance()->LASTNBDUPS<<","
            <<GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES<<","
            <<dl;
        delete geneTree;
        delete speciesTree;
        return;
    }

    Node* polygenetree = new Node(false);

    for (int i = 0; i < highspecs.size(); i++)
    {
        Node* ti = new Node(false);
        ti->CopyFrom(highspecs[i]);
        polygenetree->AddSubTree(ti);
    }

    unordered_map<Node*, Node*> polymap = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(polygenetree, speciesTree, "__", 1);


    Node* res = PolySolver::Instance()->SolvePolytomies(polygenetree, speciesTree, polymap);


    unordered_map<Node*, Node*> resmap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(res, speciesTree, "__", 1);
    int dl = GeneSpeciesTreeUtil::Instance()->GetDLScore(res, speciesTree, resmap);

    cout<<GeneSpeciesTreeUtil::Instance()->LASTNBDUPS<<","
        <<GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES<<","
        <<dl;

    delete polygenetree;
    delete res;
    delete geneTree;
    delete speciesTree;

}


string DoSGTOnHighSpecs(string gcontent, string scontent, bool preserveDupSpec)
{
    Node* speciesTree = NewickLex::ParseNewickString(scontent, true);
    Node* geneTree = NewickLex::ParseNewickString(gcontent, false);


    unordered_map<Node*, Node*> lcamap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, "__", 1);


    vector<Node*> highspecs = GeneSpeciesTreeUtil::Instance()->GetGeneTreeHighestSpeciations(geneTree, speciesTree, lcamap);

    if (highspecs.size() <= 1)
    {
        //nothing to do...
        int dl = GeneSpeciesTreeUtil::Instance()->GetDLScore(geneTree, speciesTree, lcamap);

        cout<<GeneSpeciesTreeUtil::Instance()->LASTNBDUPS<<","
            <<GeneSpeciesTreeUtil::Instance()->LASTNBLOSSES<<","
            <<dl;
        delete geneTree;
        delete speciesTree;
        return gcontent;
    }


    string gstr = "";
    //DoSGT expects a list of strings, so we'll build it for him.
    //A waste of cpu time but a gain of programmer time.
    for (int i = 0; i < highspecs.size(); i++)
    {
        gstr += NewickLex::ToNewickString(highspecs[i]);
    }

    string res = DoSuperGeneTree(gstr, scontent, preserveDupSpec, "stats");

    delete geneTree;
    delete speciesTree;

    return res;
}
