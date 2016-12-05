import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

public class main {
	/*
	 * HashMap<Integer, ArrayList<Integer>> G; HashMap<Integer,
	 * ArrayList<Boolean>> X; HashMap<Integer, ArrayList<Double>> F; // u,c
	 * HashMap<Integer, ArrayList<Double>> W;// k,c
	 * 
	 */
	HashMap<Integer, Integer> nodeIdMap; // Code id to actual id
	HashMap<Integer, Integer> nodeIdMapReverse; // Actual id to id in code.

	ArrayList<HashMap<Integer, Integer>> G; // Done
	boolean[][] X; // Done
	double[][] F;
	double[][] E;
	double[] W;

	ArrayList<HashMap<Integer, Integer>> nodeAttributeMap;
	ArrayList<HashMap<Integer, Integer>> attributeNodeMap;

	int C; // No of communities
	int numAttr;
	double lambda;
	double eta;
	int maxIter;
	double delta;
	int THREADS;
	int edges;

	double alpha;
	int V;
	
	public main(int com) {
		nodeIdMap = new HashMap<>();
		nodeIdMapReverse = new HashMap<>();
		C = com;
		numAttr = 0;
		alpha = 0.5;
		eta = 0.00005;
		lambda = 1;
		maxIter = 1000;
		THREADS = 4;
		G = new ArrayList<>();
		nodeAttributeMap = new ArrayList<>();
		attributeNodeMap = new ArrayList<>();
	}

	private void readGraph(String graphfilename, String attrfilename) throws Exception {

		BufferedReader br = new BufferedReader(new FileReader(attrfilename));
		String line;
		int nodeCounter = 0;
		ArrayList<String> attributeStringList = new ArrayList<>();
		while ((line = br.readLine()) != null) {
			attributeStringList.add(line);
		}
		V = attributeStringList.size();
		X = new boolean[V][];

		for (String l : attributeStringList) {

			String[] attributes = l.split(" ");
			if (numAttr == 0) {
				numAttr = attributes.length - 1;
			}

			int v = Integer.parseInt(attributes[0]);

			if (!nodeIdMapReverse.containsKey(v)) {
				nodeIdMap.put(nodeCounter, v);
				nodeIdMapReverse.put(v, nodeCounter);
				nodeCounter++;
			}
			int x = nodeIdMapReverse.get(v);
			X[x] = new boolean[numAttr];

			for (int i = 1; i < attributes.length; i++) {
				X[x][i - 1] = Boolean.parseBoolean(attributes[i]);
			}
		}
		br.close();
		int V = X.length;

		br = new BufferedReader(new FileReader(graphfilename));

		ArrayList<String> edgeList = new ArrayList<>();

		while ((line = br.readLine()) != null) {
			edgeList.add(line);
		}

		for (int i = 0; i < V; i++)
			G.add(new HashMap<Integer, Integer>());

		for (String l : edgeList) {
			String[] edge = l.split(" ");
			int e0 = Integer.parseInt(edge[0]);
			int e1 = Integer.parseInt(edge[1]);

			int a = nodeIdMapReverse.get(e0);
			int b = nodeIdMapReverse.get(e1);

			G.get(a).put(b, 1);
			edges++;
		}

		br.close();
		// adding back nodes that are in attrfile but not in graph
		delta = -Math.log(1 - 1.0 / G.size());
		
		edges /= 2;

		W = new double[C];
		E = new double[numAttr][C];
		F = new double[V][C];
		// Init.
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
	}

	private void initAffiliations() {
		// Compute the conductances of nodes
		HashMap<Integer, Double> conductance = new HashMap<>();
		for (int node = 0; node < G.size(); node++) {
			ArrayList<Integer> neighbors = new ArrayList<Integer>();
			neighbors.add(node);
			for (int n : G.get(node).keySet()) {
				neighbors.add(n);
			}
			int outedges = 0, inedges = 0;
			for (int v : neighbors) {
				for (int v1 : G.get(v).keySet()) {
					if (neighbors.contains(v1))
						inedges++;
					else
						outedges++;
				}
			}
			if (inedges == 0 && outedges == 0) {
				conductance.put(node, Double.MAX_VALUE);
			} else {
				conductance.put(node, (double) (outedges * 1.0 / inedges));
			}
		}

		conductance = (HashMap<Integer, Double>) mapUtil.sortByValue(conductance);
		ArrayList<Integer> alreadyAdded = new ArrayList<Integer>();
		int community = 0;
		double bias = 0.1;
		Random r = new Random();

		HashMap<Integer, ArrayList<Integer>> initialCommunities = new HashMap<>();
		// Now initialize communities based on the conductance values computed
		for (int k : conductance.keySet()) {
			// if (conductance.get(k) < 0)
			// continue;
			if (alreadyAdded.contains(k))
				continue;
			ArrayList<Integer> neighbors = new ArrayList<>();
			neighbors.add(k);
			for (int n : G.get(k).keySet()) {
				neighbors.add(n);
			}
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

		// Initialize weights E k,c, Wc

		for (int c = 0; c < C; c++) {
			double density = 0.0;
			double maxEdges = initialCommunities.get(c).size() * (initialCommunities.get(c).size() - 1);
			for (int n : initialCommunities.get(c)) {
				for (int n2 : initialCommunities.get(c)) {
					if (G.get(n).containsKey(n2)) {
						density += 1;
					}
				}
			}
			density /= maxEdges;

			W[c] = density;
		}

		for (int attr = 0; attr < numAttr; attr++) {
			for (int com = 0; com < C; com++) {
				E[attr][com] =  (r.nextDouble() * 1.0) - 0.5;
			}
		}

	}

	double computePhi(int u, int v, int c, double[] W, double[][] E) {
		double sum = 0.0;
		for (int k : nodeAttributeMap.get(u).keySet()) {
			if (nodeAttributeMap.get(v).containsKey(k))
				sum += E[k][c];
		}
		return alpha * W[c] + (1 - alpha) * sum;
	}

	double computePsi(int u, int v, double[] Fu, double[] Fv, double[] W, double[][] E) {
		double val = 0.0;
		for (int c = 0; c < C; c++) {
			double val_c = Fu[c] * Fv[c] * computePhi(u, v, c, W, E);
			val += val_c;
		}
		return val;
	}

	// L_G
	double likelihood(double[][] F, double[][] E, double[] W) {
		double Lg = 0.0;
		for (int u = 0; u < G.size(); u++) {
			for (int v = 0; v < G.size(); v++) {
				if (G.get(u).containsKey(v))
					Lg += Math.log(1 - Math.exp(-computePsi(u, v, F[u], F[v], W, E)));
				else
					Lg -= computePsi(u, v, F[u], F[v], W, E);

			}
		}
		Lg -= lambda * MathFunctions.L1Norm(E);
		return Lg;
	}

	void updateFDriver() throws InterruptedException {
		double[] sumFvc = new double[C];
		double[][] updates = new double[V][C];

		for (int v = 0; v < G.size(); v++) {
			double[] Fv = F[v];
			for (int c = 0; c < C; c++) {
				sumFvc[c] += Fv[c];
			}
		}
		int splitLen = (int) Math.ceil(G.size() * 1.0 / THREADS);
		int[] startIndices = new int[THREADS];
		int[] endIndices = new int[THREADS];

		for (int i = 0; i < THREADS; i++) {
			startIndices[i] = splitLen * i;
			endIndices[i] = Math.min(splitLen * (i + 1) - 1, V - 1);
		}

		UpdateFThread[] threads = new UpdateFThread[THREADS];

		for (int j = 0; j < threads.length; j++) {
			threads[j] = new UpdateFThread(sumFvc, startIndices[j], endIndices[j], updates);
		}

		for (UpdateFThread thread : threads) {
			thread.start();
		}
		for (UpdateFThread thread : threads) {
			thread.join();
		}

		for (int j = 0; j < THREADS; j++) {

			for (int u = startIndices[j]; u <= endIndices[j]; u++) {
				for (int c = 0; c < C; c++) {
					updates[u][c] = threads[j].updates[u][c];
				}
			}
		}
		// eta = 0.001;
		// eta = lineSearch(updates);

		for (int u = 0; u < G.size(); u++) {
			for (int c = 0; c < C; c++) {
				F[u][c] = Math.max(0.0, F[u][c] + eta * updates[u][c]);
			}
		}
	}

	class UpdateFThread extends Thread {

		// ArrayList<Integer> nodes;
		double[] sumFvc;
		// HashMap<Integer, ArrayList<Double>> updates;
		int startIndex;
		int endIndex;
		double[][] updates;

		UpdateFThread(double[] sumFvc, int startIndex, int endIndex, double[][] updates) {
			this.sumFvc = sumFvc;
			this.startIndex = startIndex;
			this.endIndex = endIndex;
			this.updates = updates;
		}

		public void run() {
			updates = updateF(sumFvc, startIndex, endIndex, updates);
		}
	}

	class UpdateEThread extends Thread {

		ArrayList<Integer> attributes;
		// HashMap<Integer, ArrayList<Double>> updates;
		double[][] updates;
		int startIndex;
		int endIndex;

		UpdateEThread(int startIndex, int endIndex, double[][] updates) {
			// this.attributes = attributes;
			this.startIndex = startIndex;
			this.endIndex = endIndex;
			this.updates = updates;
		}

		public void run() {
			updates = updateE(startIndex, endIndex, updates);
		}
	}

	class UpdateWThread extends Thread {

		// HashMap<Integer, ArrayList<Double>> updates;
		double[] updates;
		int startIndex;
		int endIndex;

		UpdateWThread(int startIndex, int endIndex, double[] updates) {
			// this.attributes = attributes;
			this.startIndex = startIndex;
			this.endIndex = endIndex;
			this.updates = updates;
		}

		public void run() {
			updates = updateW(startIndex, endIndex, updates);
		}
	}

	// Done
	double[][] updateF(double[] sumFvc, int startIndex, int endIndex, double[][] updates) {
		// HashMap<Integer, ArrayList<Double>> newF = new HashMap<>();
		// HashMap<Integer, ArrayList<Double>> updates = new HashMap<>();
		/*
		 * double[] sumFvc = new double[C];
		 * 
		 * for (int v : G.keySet()) { ArrayList<Double> Fv = F.get(v); for (int
		 * c = 0; c < C; c++) { sumFvc[c] += Fv.get(c); } }
		 */

		for (int u = startIndex; u <= endIndex; u++) {
			// newF.put(u, new ArrayList<Double>());
			// updates.put(u, new ArrayList<Double>());
			double[] derivative_G = new double[C];
			double[] derivative_X = new double[C];
			double[] diff = new double[C];
			for (int v : G.get(u).keySet()) {

				double psi = computePsi(u, v, F[u], F[v], W, E);

				double exp = Math.exp(-1 * psi);
				double fraction = exp / (1 - exp);

				for (int c = 0; c < C; c++) {
					derivative_G[c] += F[v][c] * fraction * computePhi(u, v, c, W, E);
					diff[c] -= F[v][c];
					// derivative_G[c] -= F.get(v).get(c);
				}

			}
			for (int c = 0; c < C; c++) {
				derivative_G[c] -= (diff[c] + sumFvc[c] - F[u][c]) * alpha * W[c];
				double nonOptTerm = 0.0;
				for (int k : nodeAttributeMap.get(u).keySet()) {
					for (int v : attributeNodeMap.get(k).keySet())
						if (!G.get(u).containsKey(v))
							nonOptTerm += E[k][c];
				}
				derivative_G[c] -= (1 - alpha) * nonOptTerm;
			}

			for (int c = 0; c < C; c++) {
				double update = derivative_G[c];
				updates[u][c] = update;
				// updates.get(u).add(update);
			}
		}

		// F = new HashMap<>(newF);
		// System.out.println(F);
		return updates;
	}

	double[][] updateE(int startIndex, int endIndex, double[][] updates) {

		for (int k = startIndex; k <= endIndex; k++) {
			double update = 0.0;
			for (int c = 0; c < C; c++) {
				ArrayList<Integer> nodes = new ArrayList<>(attributeNodeMap.get(k).keySet());
				for (int i = 0; i < nodes.size(); i++) {
					int u = nodes.get(i);
					for (int j = i + 1; j < nodes.size(); j++) {
						int v = nodes.get(j);
						if (G.get(u).containsKey(v)) {
							double psi = computePsi(u, v, F[u], F[v], W, E);
							double exp = Math.exp(-1 * psi);
							double fraction = exp / (1 - exp);
							update += fraction * F[u][c] * F[v][c];
						} else {
							update -= F[u][c] * F[v][c];
						}
					}

				}
				update *= (1 - alpha);
				if (E[k][c] > 0)
					update -= lambda;
				else
					update += lambda;
				updates[k][c] = update;
			}
		}
		return updates;
	}

	double[] updateW(int startIndex, int endIndex, double[] updates) {
		for (int c = startIndex; c <= endIndex; c++) {
			double update = 0.0;
			for (int u = 0; u < G.size(); u++) {
				for (int v = u + 1; v < G.size(); v++) {
					if (G.get(u).containsKey(v)) {
						double psi = computePsi(u, v, F[u], F[v], W, E);
						double exp = Math.exp(-1 * psi);
						double fraction = exp / (1 - exp);
						update += fraction * F[u][c] * F[v][c];
					} else
						update -= F[u][c] * F[v][c];

				}
			}
			update *= alpha;
			updates[c] = update;
		}
		return updates;
	}

	void updateWDriver() throws InterruptedException {
		int splitLen = (int) Math.ceil(C * 1.0 / THREADS);
		int[] startIndices = new int[THREADS];
		int[] endIndices = new int[THREADS];

		for (int i = 0; i < THREADS; i++) {
			startIndices[i] = splitLen * i;
			endIndices[i] = Math.min(splitLen * (i + 1) - 1, C - 1);
		}

		UpdateWThread[] threads = new UpdateWThread[THREADS];
		double[] updates = new double[C];

		for (int j = 0; j < threads.length; j++) {
			threads[j] = new UpdateWThread(startIndices[j], endIndices[j], updates);
		}

		for (UpdateWThread thread : threads) {
			thread.start();
		}
		for (UpdateWThread thread : threads) {
			thread.join();
		}

		// HashMap<Integer, ArrayList<Double>> updates = new HashMap<>();

		for (int j = 0; j < THREADS; j++) {
			for (int c = startIndices[j]; c <= endIndices[j]; c++) {
				updates[c] = threads[j].updates[c];
			}
		}
		// eta = lineSearchW(updates);
		// eta = 0.0001;

		for (int c = 0; c < C; c++) {
			W[c] = W[c] + eta * updates[c];
		}

	}

	void updateEDriver() throws InterruptedException {

		int splitLen = (int) Math.ceil(numAttr * 1.0 / THREADS);
		int[] startIndices = new int[THREADS];
		int[] endIndices = new int[THREADS];

		for (int i = 0; i < THREADS; i++) {
			startIndices[i] = splitLen * i;
			endIndices[i] = Math.min(splitLen * (i + 1) - 1, numAttr - 1);
		}

		UpdateEThread[] threads = new UpdateEThread[THREADS];
		double[][] updates = new double[numAttr][C];

		for (int j = 0; j < threads.length; j++) {
			threads[j] = new UpdateEThread(startIndices[j], endIndices[j], updates);
		}

		for (UpdateEThread thread : threads) {
			thread.start();
		}
		for (UpdateEThread thread : threads) {
			thread.join();
		}

		// HashMap<Integer, ArrayList<Double>> updates = new HashMap<>();

		for (int j = 0; j < THREADS; j++) {
			for (int k = startIndices[j]; k <= endIndices[j]; k++) {
				for (int c = 0; c < C; c++) {
					updates[k][c] = threads[j].updates[k][c];
				}
			}
		}
		// eta = lineSearchW(updates);
		// eta = 0.0001;

		for (int k = 0; k < numAttr; k++) {
			for (int c = 0; c < C; c++) {
				E[k][c] = E[k][c] + eta * updates[k][c];
			}
		}

	}

	void gradientAscent() throws InterruptedException {
		double startTime = System.currentTimeMillis();

		double previousLikelihood = -1;
		for (int iter = 0; iter < maxIter; iter++) {
			double likelihood = likelihood(F, E, W);
			System.out.println("Likelihood at iter " + iter + " " + likelihood);
			if (previousLikelihood != -1) {
				// System.out.println((likelihood - previousLikelihood) /
				// previousLikelihood);
				if ((previousLikelihood - likelihood) / previousLikelihood < 0.0001)
					break;
			}
			previousLikelihood = likelihood;

			updateFDriver();
			// System.out.println("Likelihood after F update at iter " + iter +
			// " " + likelihood(F, W));
			updateEDriver();
			updateWDriver();

			// System.out.println("Time " + (System.currentTimeMillis() -
			// startTime) / (1000.0 * (iter + 1)));
		}
		System.out.println("Time " + (System.currentTimeMillis() - startTime) / (1000.0));
	}

	void getCommunities() {
		ArrayList<ArrayList<Integer>> communities = new ArrayList<>();
		for (int c = 0; c < C; c++) {
			communities.add(new ArrayList<>());
		}

		for (int u = 0; u < G.size(); u++) {
			for (int c = 0; c < C; c++) {
				if (F[u][c] > delta) {
					communities.get(c).add(u);
				}
			}
		}

		for (ArrayList<Integer> comm : communities) {
			System.out.print("Circle");
			for (int x : comm)
				System.out.print(" " + nodeIdMap.get(x));
			System.out.println();
		}
	}

	int numberOfParameters() {
		return G.size() * C + numAttr * C;
	}

	double driver(String graphfilename, String attrfilename) throws Exception {
		// String graphfilename = "facebook/0.edges";
		// String attrfilename = "facebook/0.feat";
		readGraph(graphfilename, attrfilename);
		initAffiliations();
		gradientAscent();
		getCommunities();
//		for (int k = 0; k < numAttr; k++) {
//			for (int c = 0; c < C; c++) {
//				System.out.print(E[k][c] + " ");
//			}
//			System.out.println();
//		}
		return 0.0;
		// return likelihood(F, W);
	}

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		main m = new main(10);
		m.driver("facebook/107.edges", "facebook/107.feat");
	}

}
