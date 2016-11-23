import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class main {
	HashMap<Integer, ArrayList<Integer>> G;
	HashMap<Integer, ArrayList<Boolean>> X;
	HashMap<Integer, ArrayList<Double>> F; // u,c
	HashMap<Integer, ArrayList<Double>> W;// k,c
	int C; // No of communities

	public main(int C) {
		G = new HashMap<>();
		X = new HashMap<>();
		F = new HashMap<>();
		W = new HashMap<>();
	}

	private void readGraph(String graphfilename, String attrfilename) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(graphfilename));
		String line;
		while ((line = br.readLine()) != null) {
			String[] edge = line.split(" ");
			int e0 = Integer.parseInt(edge[0]);
			int e1 = Integer.parseInt(edge[1]);
			if (!G.containsKey(e0)) {
				G.put(e0, new ArrayList<Integer>());
			}
			G.get(e0).add(e1);
		}
		br = new BufferedReader(new FileReader(attrfilename));
		while ((line = br.readLine()) != null) {
			String[] attributes = line.split(" ");
			int v = Integer.parseInt(attributes[0]);
			X.put(v, new ArrayList<Boolean>());
			for (int i = 1; i < attributes.length; i++) {
				X.get(v).add(Boolean.parseBoolean(attributes[i]));
			}
		}
		// adding back nodes that are in attrfile but not in graph
		for (int node : X.keySet()) {
			if (!G.containsKey(node)) {
				G.put(node, new ArrayList<Integer>());
			}
		}
	}

	private void initAffiliations() {
		HashMap<Integer, Double> conductance = new HashMap<>();
		for (int node : G.keySet()) {
			ArrayList<Integer> neighbors = new ArrayList<Integer>();
			neighbors.add(node);
			for (int n : G.get(node)) {
				neighbors.add(n);
			}
			int outedges = 0, inedges = 0;
			for (int v : neighbors) {
				for (int v1 : G.get(v)) {
					if (neighbors.contains(v1))
						inedges++;
					else
						outedges++;
				}
			}
			if (inedges == 0 && outedges == 0) {
				conductance.put(node, (double) -1);
			} else {
				conductance.put(node, (double) (outedges * 1.0 / inedges));
				System.out.println(inedges + " " + outedges);
			}
		}
		mapUtil M = new mapUtil();
		conductance = (HashMap<Integer, Double>) M.sortByValue(conductance);
		ArrayList<Integer> alreadyAdded = new ArrayList<Integer>();
		int community = 0;
		double bias = 0.3;
		Random r = new Random();
		for (int k : conductance.keySet()) {
			if (conductance.get(k) < 0)
				continue;
			if (alreadyAdded.contains(k))
				continue;
			ArrayList<Integer> neighbors = new ArrayList<>();
			neighbors.add(k);
			for (int n : G.get(k)) {
				neighbors.add(n);
			}
			for (int n : neighbors) {
				ArrayList<Double> Fu = new ArrayList<>();
				for (int c = 0; c < C; c++) {
					if (c == community) {
						Fu.add(0.7 + (0.3 * r.nextDouble()));
					} else {
						Fu.add(0.3 * r.nextDouble());
					}
				}
				F.put(n, Fu);
				alreadyAdded.add(n);
			}
			community++;
		}
	}

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		main m = new main(10);
		String graphfilename = "facebook/0.edges";
		String attrfilename = "facebook/0.feat";
		m.readGraph(graphfilename, attrfilename);
		m.initAffiliations();
	}

}
