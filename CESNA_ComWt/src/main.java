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
	ArrayList<ArrayList<Integer>> heldOutEdges;
	ArrayList<HashMap<Integer, Integer>> G; // Done
	boolean[][] X; // Done
	double[][] F;
	double[] W;
	double[][] cosineSim;
	ArrayList<HashMap<Integer, Integer>> nodeAttributeMap;
	ArrayList<HashMap<Integer, Integer>> attributeNodeMap;

	double stoppingCondition;
	double heldOutPercent;
	int C; // No of communities
	int numAttr;
	double lambda;
	double eta;
	int maxIter;
	static int graphNum;
	int THREADS;
	int edges;
	int nRandLim;

	double alpha;
	int V;

	public main(int com) {
		C = com;
		numAttr = 0;
		alpha = 0.9;
		stoppingCondition = 0.01;
		nRandLim = 3;
		heldOutPercent = 0.2;
		lambda = 1;
		maxIter = 2000;
		THREADS = 8;
	}

	private void performHoldOut(boolean remove) {
		heldOutEdges = new ArrayList<>();
		/// HOLD OUT SET - FOR CROSS VALIDATION////

		// AT this point start removing 20% edges to the held Out set
		if (remove) {
			int numEdgesRemoved = (int) (heldOutPercent * edges);
			// System.out.println(edges);
			// System.out.println(numEdgesRemoved);
			Random r = new Random();
			for (int er = 0; er < numEdgesRemoved; er++) {
				int randV1;
				do {
					randV1 = r.nextInt(G.size());
				} while (G.get(randV1).size() == 0);
				ArrayList<Integer> neighbors = new ArrayList<>();
				for (int n : G.get(randV1).keySet()) {
					neighbors.add(n);
				}

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

		///////////////////////////////////////////////////////////////////////////////
	}

	private double heldOutLikelihood() {
		double likelihood = 0.0;
		int skippedEdges = 0;
		// System.out.println("Held size: " + heldOutEdges.size());
		for (ArrayList<Integer> e : heldOutEdges) {
			int v1 = e.get(0);
			int v2 = e.get(1);
			double l = Math.log(1 - Math.exp(-computePsi(v1, v2, F[v1], F[v2], W)));
			// System.out.println(-computePsi(v1, v2, F[v1], F[v2], W) + " " + l
			// + " DP "
			// + MathFunctions.dotProduct(F[v1], F[v2]) + " cosineSim " +
			// cosineSim[v1][v2]);
			if (MathFunctions.dotProduct(F[v1], F[v2]) == 0 && cosineSim[v1][v2] == 0) {
				skippedEdges++;
			} else
				likelihood += l;
		}
		// System.out.println("skipped " + skippedEdges);
		return likelihood * (heldOutEdges.size() / (heldOutEdges.size() - skippedEdges));
	}

	private void readGraph(String graphfilename, String attrfilename) throws Exception {
		nodeIdMap = new HashMap<>();
		nodeIdMapReverse = new HashMap<>();
		G = new ArrayList<>();
		nodeAttributeMap = new ArrayList<>();
		attributeNodeMap = new ArrayList<>();
		BufferedReader br = new BufferedReader(new FileReader(attrfilename));
		String line;
		int nodeCounter = 0;
		ArrayList<String> attributeStringList = new ArrayList<>();
		while ((line = br.readLine()) != null) {
			attributeStringList.add(line);
		}
		V = attributeStringList.size();
		X = new boolean[V][];
		cosineSim = new double[V][V];
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
				if (Integer.parseInt(attributes[i]) == 1)
					X[x][i - 1] = true;
				else
					X[x][i - 1] = false;
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
			G.get(b).put(a, 1);
		}
		edges = 0;
		for (int v = 0; v < G.size(); v++) {
			edges += G.get(v).size();
		}
		edges = edges / 2;
		br.close();

		// adding back nodes that are in attrfile but not in graph

		W = new double[C];
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
				else {
					cosineSim[u][v] = (att_v.size() * 1.0) / (Math.sqrt(uSize) * Math.sqrt(vSize));
				}

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

		// Computing the constant sim term part of likelihood

		for (int u = 0; u < G.size(); u++) {
			for (int v = 0; v < G.size(); v++) {
				if (!G.get(u).containsKey(v)) {

				}
			}
		}
	}

	double computePsi(int u, int v, double[] Fu, double[] Fv, double[] W) {
		double val = 0.0;
		for (int c = 0; c < C; c++) {
			val += Fu[c] * Fv[c] * W[c] * alpha;
		}
		val += (1.0 - alpha) * cosineSim[u][v];
		return val;
	}

	double likelihood(double[][] F, double[] W) {

		double Lg = 0.0;
		for (int u = 0; u < G.size(); u++) {
			for (int v = u + 1; v < G.size(); v++) {
				if (G.get(u).containsKey(v))
					Lg += Math.log(1 - Math.exp(-computePsi(u, v, F[u], F[v], W)));
				else
					Lg -= computePsi(u, v, F[u], F[v], W);
			}
		}
		return Lg;
	}

	double likelihoodU(int u, double[] Fu, double[][] F, double[] W) {
		double Lgu = 0.0;
		for (int v = 0; v < G.size(); v++) {
			if (v == u)
				continue;
			if (G.get(u).containsKey(v))
				Lgu += Math.log(1 - Math.exp(-computePsi(u, v, Fu, F[v], W)));
			else
				Lgu -= computePsi(u, v, Fu, F[v], W);
		}
		return Lgu;
	}

	double lineSearchFu(int u, double[] updates) {
		double beta = 0.5;
		double fx = likelihoodU(u, F[u], F, W);
		double t = 0.1;
		double controlParameter = 0.01;
		int iter = 0;
		boolean toContinue = true;
		do {
			iter++;
			if (iter > 10)
				return 0.0;
			double[] Fu_temp = new double[C];

			for (int c = 0; c < C; c++) {
				Fu_temp[c] = Math.max(0, F[u][c] + (t * updates[c]));
			}

			double fx_new = likelihoodU(u, Fu_temp, F, W);
			if (!Double.isFinite(fx_new)) {
				t = beta * t;
				continue;
			}
			if ((fx_new - fx) < t * Math.sqrt(MathFunctions.dotProduct(updates, updates)) * controlParameter) {
				// System.out.println(t);
				t = beta * t;
				continue;
			} else
				toContinue = false;
		} while (toContinue);

		return t;
	}

	double lineSearchW(double[] updates) {
		double beta = 0.5;
		double fx = likelihood(F, W);
		double t = 0.1;
		double controlParameter = 0.01;
		int iter = 0;
		boolean toContinue = true;
		do {
			t = beta * t;
			iter++;
			if (iter > 10)
				return 0.0;
			double[] W_copy = new double[C];
			for (int u = 0; u < G.size(); u++) {
				for (int c = 0; c < C; c++) {
					W_copy[c] = Math.max(0, W[c] + (t * updates[c]));
				}
			}
			double fx_new = likelihood(F, W_copy);
			if (((fx_new - fx) < t * Math.sqrt(MathFunctions.dotProduct(updates, updates)) * controlParameter)
					|| !Double.isFinite(fx_new)) {
				continue;
			} else {
				// System.out.println("eta is " + t + " likelihood " +
				// likelihood(F, W_copy));
				break;
			}
		} while (toContinue);
		// System.out.println("t returned " + t);
		return t;
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

		for (int u = 0; u < G.size(); u++) {
			eta = lineSearchFu(u, updates[u]);
			for (int c = 0; c < C; c++) {
				F[u][c] = Math.max(0.0, F[u][c] + eta * updates[u][c]);
			}
		}
		// System.out.println("after update with eta = " + eta + " " +
		// likelihood(F, W));
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
			double[] diff = new double[C];
			for (int v : G.get(u).keySet()) {

				double psi = computePsi(u, v, F[u], F[v], W);

				double exp = Math.exp(-1 * psi);
				double fraction = exp / (1 - exp);

				for (int c = 0; c < C; c++) {
					derivative_G[c] += F[v][c] * fraction * alpha * W[c];
					diff[c] += F[v][c];
					// derivative_G[c] -= F.get(v).get(c);
				}

			}
			for (int c = 0; c < C; c++) {
				derivative_G[c] -= (sumFvc[c] - diff[c] - F[u][c]) * alpha * W[c];
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

	double[] updateW(int startIndex, int endIndex, double[] updates) {
		for (int c = startIndex; c <= endIndex; c++) {
			double update = 0.0;
			for (int u = 0; u < G.size(); u++) {
				for (int v = u + 1; v < G.size(); v++) {
					if (G.get(u).containsKey(v)) {
						double psi = computePsi(u, v, F[u], F[v], W);
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
		eta = lineSearchW(updates);
		// eta = 0.0001;
		// System.out.println("eta chosen for W " + eta);
		for (int c = 0; c < C; c++) {
			W[c] = Math.max(0.0, W[c] + eta * updates[c]);
		}
		// System.out.println("Likelihood after W update " + likelihood(F, W));
	}

	void gradientAscent() throws InterruptedException {
		double startTime = System.currentTimeMillis();

		double previousLikelihood = -1;
		for (int iter = 0; iter < maxIter; iter++) {
			double likelihood = likelihood(F, W);
			// System.out.println("Likelihood at iter " + iter + " " +
			// likelihood);
			if (previousLikelihood != -1) {
				// System.out.println((likelihood - previousLikelihood) /
				// previousLikelihood);
				if ((previousLikelihood - likelihood) / previousLikelihood < stoppingCondition)
					break;
			}
			previousLikelihood = likelihood;

			updateFDriver();
			// System.out.println("Likelihood after F update at iter " + iter +
			// " " + likelihood(F, W));
			// updateEDriver();
			// updateWDriver();

			// System.out.println("Time " + (System.currentTimeMillis() -
			// startTime) / (1000.0 * (iter + 1)));
		}
		// System.out.println("Time " + (System.currentTimeMillis() - startTime)
		// / (1000.0));
	}

	void getCommunities(boolean print) {
		double[] delta_c = new double[C];
		ArrayList<ArrayList<Integer>> communities = new ArrayList<>();
		for (int c = 0; c < C; c++) {
			communities.add(new ArrayList<>());
			// System.out.println("w_c " + W[c]);
			delta_c[c] = Math.sqrt(Math.log(G.size() * 1.0 / (G.size() - 1)) / (W[c] * alpha * 1.0));
			// System.out.println(Math.log((G.size() * 1.0) / (G.size() - 1)));
			// System.out.println("hello " + Math.log(G.size() * 1.0 / (G.size()
			// - 1)) / (W[c] * alpha * 1.0));
			// System.out.println("delta_c " + delta_c[c]);
		}
		// System.out.println(delta_c);

		for (int u = 0; u < G.size(); u++) {
			for (int c = 0; c < C; c++) {
				if (F[u][c] > delta_c[c]) {
					// System.out.println(F[u][c]);
					// System.out.println(delta_c[c]);
					communities.get(c).add(u);
				}
			}
		}
		if (print) {
			for (ArrayList<Integer> comm : communities) {
				System.out.print("Circle");
				for (int x : comm)
					System.out.print(" " + nodeIdMap.get(x));
				System.out.println();
			}
		}
	}

	int numberOfParameters() {
		return G.size() * C + numAttr * C;
	}

	double driver(String graphFileName, String attrFileName, boolean print) throws Exception {
		// String graphfilename = "facebook/0.edges";
		// String attrfilename = "facebook/0.feat";

		initAffiliations();
		gradientAscent();
		getCommunities(print);
		// for (int k = 0; k < numAttr; k++) {
		// for (int c = 0; c < C; c++) {
		// System.out.print(E[k][c] + " ");
		// }
		// System.out.println();
		// }
		return 0.0;
		// return likelihood(F, W);
	}

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		main m = new main(0);
		graphNum = 3980;
		HashMap<Integer, Double> perf = new HashMap<>();
		String graphFileName = "facebook/" + String.valueOf(graphNum) + ".edges";
		String attrFileName = "facebook/" + String.valueOf(graphNum) + ".feat";
		int bestK = 0;
		double bestKLikelihood = (double) 0;
		m.readGraph(graphFileName, attrFileName);
		System.out.println(m.V);
		int lowerLimit = 4;
		int upperLimit = 8;
		if (m.V >= 400){
			lowerLimit = 6;
			upperLimit = 14;
		}
		for (int K = lowerLimit; K <= upperLimit; K += 2) {
			System.out.println(K);
			m.C = K;
			double avgL = 0.0;
			for (int nRand = 0; nRand < m.nRandLim; nRand++) {
				System.out.println(nRand);
				m.readGraph(graphFileName, attrFileName);
				m.performHoldOut(true);
				m.driver(graphFileName, attrFileName, false);
				double L = m.heldOutLikelihood();
				avgL += L;
				m.performHoldOut(false);
			}
			avgL /= m.nRandLim;
			perf.put(K, avgL);
			if (bestKLikelihood == 0.0) {
				bestKLikelihood = avgL;
				bestK = K;
			} else {
				if (avgL > bestKLikelihood) {
					bestKLikelihood = avgL;
					bestK = K;
				}
			}
		}
		m.C = bestK;
		m.stoppingCondition /= 1000;
		for (int i = 0; i < 5; i++){
			m.readGraph(graphFileName, attrFileName);
			m.driver(graphFileName, attrFileName, true);
			System.out.println("-");
		}
		System.out.println(bestK);
		System.out.println(perf);
	}

}
