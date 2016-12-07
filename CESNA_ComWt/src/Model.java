import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Model {
	/*
	 * Mapping node ids from input graph to the range [0,V).
	 */
	HashMap<Integer, Integer> nodeIdMap; // Code id to actual id
	HashMap<Integer, Integer> nodeIdMapReverse; // Actual id to id in code
	static ArrayList<HashMap<Integer, Integer>> G; // Graph stored as adj. list
	boolean[][] X; // Boolean matrix to store node attributes
	static double[][] cosineSim; // Pre-computed attribute cosine similarity
	static double alpha; // Hyper-parameter

	double[][] F; // User community memberships
	double[] W; // Community density
	// Held-out edges for choosing communities.
	ArrayList<ArrayList<Integer>> heldOutEdges;
	// Attributes present for a node.
	ArrayList<HashMap<Integer, Integer>> nodeAttributeMap;
	// Nodes containing that attribute.
	ArrayList<HashMap<Integer, Integer>> attributeNodeMap;

	// Relative change in log likelihood for convergence
	double stoppingCondition;
	// Percentage of edges to remove for choosing communities.
	double heldOutPercent;
	int C; // # of communities
	int numAttr; // # of attributes for each node
	double eta; // Learning rate (set by line search)
	int maxIter; // Maximum iterations for the gradient ascent algorithm
	int THREADS; // Threads for parallelism
	int V; // # Nodes in the graph
	int edges; // # Edges in the graph

	public Model(int com) {
		// Initialize basic parameters for the model.
		C = com;
		numAttr = 0; // Will be set later
		alpha = 0.9;
		stoppingCondition = 0.01;
		heldOutPercent = 0.2;
		maxIter = 2000;
		THREADS = 8;
	}

	/**
	 * Function to hold out a set of edges to choose the best number of
	 * communities.
	 */
	void performHoldOut(boolean remove) {
		heldOutEdges = new ArrayList<>();
		if (remove) {
			int numEdgesRemoved = (int) (heldOutPercent * edges);
			Random r = new Random();
			for (int er = 0; er < numEdgesRemoved; er++) {
				int randV1;
				do {
					randV1 = r.nextInt(G.size());
				} while (G.get(randV1).size() == 0);
				ArrayList<Integer> neighbors = new ArrayList<>();
				for (int n : G.get(randV1).keySet())
					neighbors.add(n);
				int randV2;
				randV2 = r.nextInt(G.get(randV1).size());
				ArrayList<Integer> e = new ArrayList<>();
				e.add(randV1);
				int v2 = neighbors.get(randV2);
				e.add(v2);
				heldOutEdges.add(e);
				G.get(randV1).remove(v2);
				G.get(v2).remove(randV1);
				edges--;
			}
		} else {
			for (ArrayList<Integer> e : heldOutEdges) {
				edges++;
				G.get(e.get(0)).put(e.get(1), 0);
				G.get(e.get(1)).put(e.get(0), 0);
			}
			heldOutEdges = new ArrayList<>();
		}
	}

	/**
	 * Function to read the input graph and node attributes, and then initialize
	 * the pre-computed model parameters - The cosine similarity b/w every pair
	 * of users in pre-computed here.
	 */
	void readGraph(String graphfilename, String attrfilename) throws Exception {
		nodeIdMap = new HashMap<>();
		nodeIdMapReverse = new HashMap<>();
		G = new ArrayList<>();
		nodeAttributeMap = new ArrayList<>();
		attributeNodeMap = new ArrayList<>();
		// Read the node attributes first. Here, we map each node id to the
		// range [0,V)
		BufferedReader br = new BufferedReader(new FileReader(attrfilename));
		String line;
		int nodeCounter = 0;
		ArrayList<String> attributeStringList = new ArrayList<>();
		while ((line = br.readLine()) != null)
			attributeStringList.add(line);
		V = attributeStringList.size();
		X = new boolean[V][];
		cosineSim = new double[V][V];
		for (String l : attributeStringList) {
			String[] attributes = l.split(" ");
			if (numAttr == 0)
				numAttr = attributes.length - 1;
			int v = Integer.parseInt(attributes[0]);
			if (!nodeIdMapReverse.containsKey(v)) {
				nodeIdMap.put(nodeCounter, v);
				nodeIdMapReverse.put(v, nodeCounter);
				nodeCounter++;
			}
			int x = nodeIdMapReverse.get(v);
			X[x] = new boolean[numAttr];
			for (int i = 1; i < attributes.length; i++) {
				if (Integer.parseInt(attributes[i]) == 1)
					X[x][i - 1] = true;
				else
					X[x][i - 1] = false;
			}
		}
		br.close();
		int V = X.length;
		// Read the graph next.
		br = new BufferedReader(new FileReader(graphfilename));
		ArrayList<String> edgeList = new ArrayList<>();
		while ((line = br.readLine()) != null)
			edgeList.add(line);

		for (int i = 0; i < V; i++)
			G.add(new HashMap<Integer, Integer>());

		for (String l : edgeList) {
			String[] edge = l.split(" ");
			int e0 = Integer.parseInt(edge[0]);
			int e1 = Integer.parseInt(edge[1]);
			int a = nodeIdMapReverse.get(e0);
			int b = nodeIdMapReverse.get(e1);
			G.get(a).put(b, 1);
			G.get(b).put(a, 1);
		}
		edges = 0;
		for (int v = 0; v < G.size(); v++)
			edges += G.get(v).size();
		edges = edges / 2;
		br.close();
		// Initialize W and F
		W = new double[C];
		F = new double[V][C];

		for (int k = 0; k < numAttr; k++)
			attributeNodeMap.add(new HashMap<>());

		for (int u = 0; u < V; u++) {
			nodeAttributeMap.add(new HashMap<>());
			for (int k = 0; k < numAttr; k++)
				if (X[u][k]) {
					nodeAttributeMap.get(u).put(k, 0);
					attributeNodeMap.get(k).put(u, 0);
				}
		}
		// Compute cosine similarity based on node attributes.
		for (int u = 0; u < G.size(); u++) {
			Set<Integer> att_u = nodeAttributeMap.get(u).keySet();
			int uSize = nodeAttributeMap.get(u).size();
			for (int v = 0; v < G.size(); v++) {
				HashSet<Integer> att_v = new HashSet<>();
				Set<Integer> att_v1 = nodeAttributeMap.get(v).keySet();
				for (int v1 : att_v1)
					att_v.add(v1);
				int vSize = nodeAttributeMap.get(v).size();
				att_v.retainAll(att_u);
				if (uSize == 0 || vSize == 0)
					cosineSim[u][v] = 0.0;
				else
					cosineSim[u][v] = (att_v.size() * 1.0) / (Math.sqrt(uSize) * Math.sqrt(vSize));
			}
		}
	}

	/**
	 * Initial community memberships for each (u,c) pair are assigned by
	 * computing conductances and in-out edge ratios.
	 */
	private void initAffiliations() {
		// Compute the conductances of nodes
		HashMap<Integer, Double> conductance = new HashMap<>();
		for (int node = 0; node < G.size(); node++) {
			ArrayList<Integer> neighbors = new ArrayList<Integer>();
			neighbors.add(node);
			for (int n : G.get(node).keySet())
				neighbors.add(n);

			int outedges = 0, inedges = 0;
			for (int v : neighbors) {
				for (int v1 : G.get(v).keySet()) {
					if (neighbors.contains(v1))
						inedges++;
					else
						outedges++;
				}
			}
			if (inedges == 0 && outedges == 0)
				conductance.put(node, Double.MAX_VALUE);
			else
				conductance.put(node, (double) (outedges * 1.0 / inedges));
		}

		conductance = (HashMap<Integer, Double>) mapUtil.sortByValue(conductance);
		ArrayList<Integer> alreadyAdded = new ArrayList<Integer>();
		int community = 0;
		double bias = 0.1;
		Random r = new Random();
		HashMap<Integer, ArrayList<Integer>> initialCommunities = new HashMap<>();
		// Now initialize communities based on the conductance values computed
		for (int k : conductance.keySet()) {
			if (alreadyAdded.contains(k))
				continue;
			ArrayList<Integer> neighbors = new ArrayList<>();
			neighbors.add(k);
			for (int n : G.get(k).keySet())
				neighbors.add(n);

			initialCommunities.put(community, neighbors);
			for (int n : neighbors) {
				for (int c = 0; c < C; c++) {
					if (c == community)
						F[n][c] = 1 - bias + (bias * r.nextDouble());
					else
						F[n][c] = bias * r.nextDouble();
				}
				alreadyAdded.add(n);
			}
			community++;
		}
		// Initialize W based on community density.
		for (int c = 0; c < C; c++) {
			double density = 0.0;
			double maxEdges = initialCommunities.get(c).size() * (initialCommunities.get(c).size() - 1);
			for (int n : initialCommunities.get(c)) {
				for (int n2 : initialCommunities.get(c)) {
					if (G.get(n).containsKey(n2))
						density += 1;
				}
			}
			density /= maxEdges;
			W[c] = density;
		}
	}

	/**
	 * Performs gradient update for F by creating threads for computing the
	 * update values in parallel for different users.
	 */
	void updateFDriver() throws InterruptedException {
		double[] sumFvc = new double[C];
		double[][] updates = new double[V][C];

		for (int v = 0; v < G.size(); v++) {
			double[] Fv = F[v];
			for (int c = 0; c < C; c++)
				sumFvc[c] += Fv[c];
		}
		int splitLen = (int) Math.ceil(G.size() * 1.0 / THREADS);
		int[] startIndices = new int[THREADS];
		int[] endIndices = new int[THREADS];

		for (int i = 0; i < THREADS; i++) {
			startIndices[i] = splitLen * i;
			endIndices[i] = Math.min(splitLen * (i + 1) - 1, V - 1);
		}
		UpdateFThread[] threads = new UpdateFThread[THREADS];
		for (int j = 0; j < threads.length; j++)
			threads[j] = new UpdateFThread(sumFvc, startIndices[j], endIndices[j], updates, F, W);

		for (UpdateFThread thread : threads)
			thread.start();

		for (UpdateFThread thread : threads)
			thread.join();

		for (int j = 0; j < THREADS; j++)
			for (int u = startIndices[j]; u <= endIndices[j]; u++)
				for (int c = 0; c < C; c++)
					updates[u][c] = threads[j].updates[u][c];

		for (int u = 0; u < G.size(); u++) {
			eta = LineSearch.getLearningRateFU(u, updates[u], F, W);
			for (int c = 0; c < C; c++)
				F[u][c] = Math.max(0.0, F[u][c] + eta * updates[u][c]);
		}
	}

	/**
	 * Performs gradient update for W by creating threads for computing the
	 * update values in parallel for different communities.
	 */
	void updateWDriver() throws InterruptedException {
		// Split the set of communities equally across threads.
		int splitLen = (int) Math.ceil(C * 1.0 / THREADS);
		int[] startIndices = new int[THREADS];
		int[] endIndices = new int[THREADS];

		for (int i = 0; i < THREADS; i++) {
			startIndices[i] = splitLen * i;
			endIndices[i] = Math.min(splitLen * (i + 1) - 1, C - 1);
		}

		UpdateWThread[] threads = new UpdateWThread[THREADS];
		double[] updates = new double[C];

		for (int j = 0; j < threads.length; j++)
			threads[j] = new UpdateWThread(startIndices[j], endIndices[j], updates, F, W);

		for (UpdateWThread thread : threads)
			thread.start();

		for (UpdateWThread thread : threads)
			thread.join();

		for (int j = 0; j < THREADS; j++)
			for (int c = startIndices[j]; c <= endIndices[j]; c++)
				updates[c] = threads[j].updates[c];

		eta = LineSearch.getLearningRateWC(updates, F, W);
		for (int c = 0; c < C; c++)
			W[c] = Math.max(0.0, W[c] + eta * updates[c]);
	}

	/**
	 * Gradient ascent function. Iteratively updates F (and W).
	 */
	void gradientAscent() throws InterruptedException {
		double previousLikelihood = -1;
		for (int iter = 0; iter < maxIter; iter++) {
			double likelihood = MathFunctions.likelihood(F, W);
			if (previousLikelihood != -1) {
				if ((previousLikelihood - likelihood) / previousLikelihood < stoppingCondition)
					break;
			}
			previousLikelihood = likelihood;
			updateFDriver();
			// updateWDriver();
		}
	}

	/**
	 * Post-processing task of getting hard community assignments from community
	 * memberships after the gradient ascent algorithm converges.
	 * 
	 * @throws IOException
	 */
	void getCommunities(boolean print, String outputFilename) throws IOException {
		// Compute threshold for community assignment first.
		double[] delta_c = new double[C];
		ArrayList<ArrayList<Integer>> communities = new ArrayList<>();
		for (int c = 0; c < C; c++) {
			communities.add(new ArrayList<>());
			delta_c[c] = Math.sqrt(Math.log(G.size() * 1.0 / (G.size() - 1)) / (W[c] * alpha * 1.0));
		}

		for (int u = 0; u < G.size(); u++)
			for (int c = 0; c < C; c++)
				if (F[u][c] > delta_c[c])
					communities.get(c).add(u);

		if (print) {
			for (ArrayList<Integer> comm : communities) {
				System.out.print("Circle");
				for (int x : comm)
					System.out.print(" " + nodeIdMap.get(x));
				System.out.println();
			}
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFilename));
		for (int i = 0; i < communities.size(); i++) {
			ArrayList<Integer> comm = communities.get(i);
			bw.write("Circle"+i);
			for (int x : comm)
				bw.write("\t" + nodeIdMap.get(x));
			bw.write("\n");
		}
		bw.close();
	}

	/**
	 * Driver function for the program which runs the model once for the current
	 * parameter settings.
	 */
	double driver(String graphFileName, String attrFileName, String outputFilename, boolean print) throws Exception {
		initAffiliations();
		gradientAscent();
		getCommunities(print, outputFilename);
		return 0.0;
	}

}
