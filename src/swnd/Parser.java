package swnd;

import edu.asu.emit.algorithm.graph.Graph;
import edu.asu.emit.algorithm.graph.Path;
//import edu.asu.emit.algorithm.graph.shortestpaths.DijkstraShortestPathAlg;
import edu.asu.emit.algorithm.graph.shortestpaths.YenTopKShortestPathsAlg;

import java.util.List;


public class Parser {

    private Graph graph;

    Parser(int type){
        graph = new Graph("data/v"+type);
    }

    Integer[] compute(int a, int b, int c) {

        YenTopKShortestPathsAlg yenAlg = new YenTopKShortestPathsAlg(graph);
        List<Path> shortest_paths_list = yenAlg.getShortestPaths(graph.getVertex(a), graph.getVertex(b), c);
        //System.out.println(":" + shortest_paths_list);

        //Object x = shortest_paths_list.toArray();
        //System.out.println(yenAlg.getResultList().size());

        if (shortest_paths_list.size() < c) {
            throw new ExceptionInInitializerError();
        }


        String path = shortest_paths_list.get(c - 1).getVertexList().toString();

        System.out.println("From: "+ a + " To: "+ b+ "  - " + path);

        String[] tokens = path.split(":");
        String[] items = tokens[0].replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(" ", "").split(",");

        Integer[] results = new Integer[items.length];

        for (int i = 0; i < items.length; i++) {
            try {
                results[i] = Integer.parseInt(items[i]);
                //System.out.print(results[i].toString());
                //System.out.println();
            } catch (NumberFormatException nfe) {
                System.err.println(nfe.toString());
            }
        }
        return results;


    }
}


