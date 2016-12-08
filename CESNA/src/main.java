import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class main {

	HashMap<Integer, Integer> nodeIdMap; // Code id to actual id
	HashMap<Integer, Integer> nodeIdMapReverse; // Actual id to id in code.
	static ArrayList<HashMap<Integer, Integer>> G;
	static boolean[][] X;
	double[][] F;
	double[][] W;

	int C; // No of communities
	int numAttr;
	static double alpha;
	double eta;
	static double lambda;
	int maxIter;
	double delta;
	int THREADS;
	int V;
	int E;

	public main(int com) {
		nodeIdMap = new HashMap<>();
		nodeIdMapReverse = new HashMap<>();
		C = com;
		numAttr = 0;
		alpha = 0.5;
		eta = 0.0015;
		lambda = 1;
		maxIter = 2000;
		THREADS = 8;
		G = new ArrayList<>();
	}

	/**
	 * Function to read the input graph and node attributes, and then initialize
	 * the pre-computed model parameters.
	 */
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
			G.get(a).put(b, 0);
			G.get(b).put(a, 0);
		}
		br.close();
		for (int u = 0; u < G.size(); u++)
			E += G.get(u).size();
		// adding back nodes that are in attrfile but not in graph
		delta = -Math.log(1 - 1.0 / G.size());
		E /= 2;
		W = new double[numAttr][C];
		F = new double[V][C];
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

		// Now initialize communities based on the conductance values computed
		for (int k : conductance.keySet()) {
			if (alreadyAdded.contains(k))
				continue;
			ArrayList<Integer> neighbors = new ArrayList<>();
			neighbors.add(k);
			for (int n : G.get(k).keySet()) {
				neighbors.add(n);
			}
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

		// Initialize weights W k,c
		for (int attr = 0; attr < numAttr; attr++) {
			for (int com = 0; com < C; com++) {
				W[attr][com] = (r.nextDouble() * 2.0) - 1.0;
			}
		}
	}

	// L_G - graph log likelihood
	double graphLikelihood(double[][] F) {
		double Lg = 0.0;
		for (int u = 0; u < G.size(); u++) {
			for (int v = u + 1; v < G.size(); v++) {
				if (G.get(u).keySet().contains(v))
					Lg += Math.log(1 - Math.exp(-1 * MathFunctions.dotProduct(F[u], F[v])));
				else
					Lg -= MathFunctions.dotProduct(F[u], F[v]);
			}
		}
		return Lg;
	}

	// Compute L_G - graph log likelihood for specific node u
	double graphLikelihoodU(int u, double[] Fu, double[][] F) {
		double Lgu = 0.0;
		for (int v = 0; v < G.size(); v++) {
			if (G.get(u).keySet().contains(v))
				Lgu += Math.log(1 - Math.exp(-1 * MathFunctions.dotProduct(F[u], F[v])));
			else
				Lgu -= MathFunctions.dotProduct(F[u], F[v]);
		}
		return Lgu;
	}

	// Lx
	double attributeLikelihood(double[][] W, double[][] F) {
		double Lx = 0.0;
		for (int u = 0; u < G.size(); u++) {
			for (int k = 0; k < numAttr; k++) {
				boolean x = X[u][k];
				double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(W[k], F[u])));
				int z = x ? 1 : 0;
				Lx += z * Math.log(Quk) + (1 - z) * Math.log(1 - Quk);
				if (Double.isNaN(Lx)) {
					System.out.println("Quk = " + Quk);
					System.out.println(MathFunctions.dotProduct(W[k], F[u]));
					System.exit(0);
				}
			}
		}
		return Lx;
	}

	// Compute L_X - attribute log likelihood
	double attributeLikelihoodK(int k, double[] Wk, double[][] F) {
		double Lxk = 0.0;
		for (int u = 0; u < G.size(); u++) {
			boolean x = X[u][k];
			double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(Wk, F[u])));
			int z = x ? 1 : 0;
			Lxk += z * Math.log(Quk) + (1 - z) * Math.log(1 - Quk);
			if (Double.isNaN(Lxk)) {
				System.out.println("Quk = " + Quk);
				System.out.println(MathFunctions.dotProduct(Wk, F[u]));
				System.exit(0);
			}
		}
		return Lxk;
	}

	// Compute log likelihood
	double likelihood(double[][] F, double[][] W) {
		return (1 - alpha) * graphLikelihood(F) + alpha * attributeLikelihood(W, F) - lambda * MathFunctions.L1Norm(W);
	}

	double lineSearchW(int k, double[] updates) {
		double fx = attributeLikelihoodK(k, W[k], F);
		double t = 0.1;
		double beta = 0.3;
		int iter = 0;
		do {
			iter++;
			if (iter > 10)
				return 0.0;

			double[] Wk = new double[C];

			for (int c = 0; c < C; c++) {
				Wk[c] = W[k][c] + (t * updates[c]);
			}

			double fx_new = attributeLikelihoodK(k, Wk, F);
			if (((fx_new - fx) < t * Math.sqrt(MathFunctions.dotProduct(updates, updates)) * 0.05)
					|| !Double.isFinite(fx_new)) {
				t = beta * t;
			} else
				return t;

		} while (true);
	}

	double lineSearchF(int u, double[] updates) {
		double beta = 0.3;
		double fx = graphLikelihoodU(u, F[u], F);
		double t = 0.1;
		int iter = 0;

		do {
			iter++;
			if (iter > 10)
				return 0.0;
			double[] Fu = new double[C];

			for (int c = 0; c < C; c++) {
				Fu[c] = Math.max(0, F[u][c] + (t * updates[c]));
			}

			double fx_new = graphLikelihoodU(u, Fu, F);
			if (((fx_new - fx) < t * Math.sqrt(MathFunctions.dotProduct(updates, updates)) * 0.05)
					|| !Double.isFinite(fx_new)) {
				t = beta * t;
			} else
				return t;
		} while (true);
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
			// eta = lineSearchF(u, updates[u]);
			for (int c = 0; c < C; c++)
				F[u][c] = Math.max(0.0, F[u][c] + eta * updates[u][c]);
		}
	}

	/**
	 * Performs gradient update for F by creating threads for computing the
	 * update values in parallel for different attributes.
	 */
	void updateWDriver() throws InterruptedException {
		int splitLen = (int) Math.ceil(numAttr * 1.0 / THREADS);
		int[] startIndices = new int[THREADS];
		int[] endIndices = new int[THREADS];

		for (int i = 0; i < THREADS; i++) {
			startIndices[i] = splitLen * i;
			endIndices[i] = Math.min(splitLen * (i + 1) - 1, numAttr - 1);
		}

		UpdateWThread[] threads = new UpdateWThread[THREADS];
		double[][] updates = new double[numAttr][C];

		for (int j = 0; j < threads.length; j++)
			threads[j] = new UpdateWThread(startIndices[j], endIndices[j], updates, F, W);

		for (UpdateWThread thread : threads)
			thread.start();

		for (UpdateWThread thread : threads)
			thread.join();

		// HashMap<Integer, ArrayList<Double>> updates = new HashMap<>();

		for (int j = 0; j < THREADS; j++)
			for (int k = startIndices[j]; k <= endIndices[j]; k++)
				for (int c = 0; c < C; c++)
					updates[k][c] = threads[j].updates[k][c];

		for (int k = 0; k < numAttr; k++) {
			// eta = lineSearchW(k, updates[k]);
			for (int c = 0; c < C; c++)
				W[k][c] = W[k][c] + eta * updates[k][c];
		}

	}

	void gradientAscent() throws InterruptedException {
		double startTime = System.currentTimeMillis();
		double previousLikelihood = -1;
		for (int iter = 0; iter < maxIter; iter++) {
			double likelihood = likelihood(F, W);
			System.out.println("Likelihood at iter " + iter + " " + likelihood);
			if (previousLikelihood != -1) {
				if ((previousLikelihood - likelihood) / previousLikelihood < 0.0001)
					break;
			}
			previousLikelihood = likelihood;
			updateFDriver();
			updateWDriver();

			// System.out.println("Time " + (System.currentTimeMillis() -
			// startTime) / (1000.0 * (iter + 1)));
		}
		System.out.println("Time " + (System.currentTimeMillis() - startTime) / (1000.0));
	}

	// Find hard community assignments based on threshold defined in the paper.
	void getCommunities(String outputfilename) throws IOException {
		ArrayList<ArrayList<Integer>> communities = new ArrayList<>();
		for (int c = 0; c < C; c++) 
			communities.add(new ArrayList<>());
		
		for (int u = 0; u < G.size(); u++)
			for (int c = 0; c < C; c++)
				if (F[u][c] > delta)
					communities.get(c).add(u);

		for (ArrayList<Integer> comm : communities) {
			System.out.print("Circle");
			for (int x : comm)
				System.out.print(" " + nodeIdMap.get(x));
			System.out.println();
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputfilename));
		for (int i = 0; i < communities.size(); i++) {
			ArrayList<Integer> comm = communities.get(i);
			System.out.print("Circle");
			for (int x : comm)
				System.out.print(" " + nodeIdMap.get(x));
			System.out.println();

			bw.write("Circle" + i);
			for (int x : comm)
				bw.write("\t" + nodeIdMap.get(x));
			bw.write("\n");
		}
		bw.close();
	}

	int numberOfParameters() {
		return G.size() * C + numAttr * C;
	}

	double driver(String graphfilename, String attrfilename, String outputfilename) throws Exception {
		readGraph(graphfilename, attrfilename);
		initAffiliations();
		gradientAscent();
		getCommunities(outputfilename);
		return 0.0;
	}

	public static void main(String[] args) throws Exception {
		main m = new main(8);
		m.driver("facebook/3980.edges", "facebook/3980.feat", "circles.txt");
	}
}
