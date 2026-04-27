package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.suffixarray.SuffixArray;

import java.util.ArrayList;
import java.util.HashMap;

public class BacktrackTable<T extends Matter> extends HashMap<T, BacktrackPointer> {
    private static final long serialVersionUID = 1L;
    DeNovoGraph<T> graph;

    public BacktrackTable(DeNovoGraph<T> graph) {
        this.graph = graph;
    }

    public void getReconstructions(T curNode, int score, String prefix, ArrayList<String> reconstructions) {
        getReconstructions(curNode, score, prefix, reconstructions, null);
    }

    public void getReconstructions(T curNode, int score, String prefix, ArrayList<String> reconstructions, SuffixArray sa) {
        if (sa != null && sa.search(prefix) < 0)
            return;

        BacktrackPointer pointer = this.get(curNode);
        if (pointer == null)
            return;
        if (score >= pointer.getMaxScore())
            return;
        assert (pointer != null);
        if (curNode.equals(graph.getSource()))    // source
        {
            reconstructions.add(prefix);
            return;
        }

        for (DeNovoGraph.Edge<T> edge : graph.getEdges(curNode)) {
            int edgeIndex = edge.getEdgeIndex();
            if (pointer.isSet(score, edgeIndex))
                getReconstructions(edge.getPrevNode(), score - (edge.getEdgeScore() + pointer.getNodeScore()), prefix + graph.getAASet().getAminoAcid(edgeIndex).getResidueStr(), reconstructions, sa);
        }
    }

    public String getOneReconstruction(T curNode, int score, String prefix) {
        BacktrackPointer pointer = this.get(curNode);
        if (pointer == null)
            return null;
        if (score >= pointer.getMaxScore())
            return null;
        assert (pointer != null);
        if (curNode.equals(graph.getSource()))    // source
        {
            return prefix;
        }
        for (DeNovoGraph.Edge<T> edge : graph.getEdges(curNode)) {
            int edgeIndex = edge.getEdgeIndex();
            if (pointer.isSet(score, edgeIndex))
                getOneReconstruction(edge.getPrevNode(), score - (edge.getEdgeScore() + pointer.getNodeScore()), prefix + graph.getAASet().getAminoAcid(edgeIndex).getResidueStr());
        }
        return null;
    }
}
