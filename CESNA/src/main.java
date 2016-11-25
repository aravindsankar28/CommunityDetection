import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;

public class main {
	/*
	 * HashMap<Integer, ArrayList<Integer>> G; HashMap<Integer,
	 * ArrayList<Boolean>> X; HashMap<Integer, ArrayList<Double>> F; // u,c
	 * HashMap<Integer, ArrayList<Double>> W;// k,c
	 * 
	 */
	HashMap<Integer, Integer> nodeIdMap; // Code id to actual id
	HashMap<Integer, Integer> nodeIdMapReverse; // Actual id to id in code.

	ArrayList<ArrayList<Integer>> G; // Done
	boolean[][] X; // Done
	double[][] F;
	double[][] W;

	int C; // No of communities
	int numAttr;
	double alpha;
	double eta;
	double lambda;
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
		eta = 0.001;
		lambda = 1;
		maxIter = 1000;
		THREADS = 4;
		G = new ArrayList<>();
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
			G.add(new ArrayList<Integer>());

		for (String l : edgeList) {
			String[] edge = l.split(" ");
			int e0 = Integer.parseInt(edge[0]);
			int e1 = Integer.parseInt(edge[1]);

			int a = nodeIdMapReverse.get(e0);
			int b = nodeIdMapReverse.get(e1);

			G.get(a).add(b);
			E++;
		}

		br.close();
		// adding back nodes that are in attrfile but not in graph
		delta = -Math.log(1 - 1.0 / G.size());
		E /= 2;

		W = new double[numAttr][C];
		F = new double[V][C];
	}

	private void initAffiliations() {
		// Compute the conductances of nodes
		HashMap<Integer, Double> conductance = new HashMap<>();
		for (int node = 0; node < G.size(); node++) {
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

		// Now initialize communities based on the conductance values computed
		for (int k : conductance.keySet()) {
			// if (conductance.get(k) < 0)
			// continue;
			if (alreadyAdded.contains(k))
				continue;
			ArrayList<Integer> neighbors = new ArrayList<>();
			neighbors.add(k);
			for (int n : G.get(k)) {
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

	// L_G
	double graphLikelihood(double[][] F) {
		double Lg = 0.0;
		for (int u = 0; u < G.size(); u++) {
			for (int v : G.get(u)) {
				if (G.get(u).contains(v))
					Lg += Math.log(1 - Math.exp(-1 * MathFunctions.dotProduct(F[u], F[v])));
				else
					Lg -= MathFunctions.dotProduct(F[u], F[v]);
			}
		}
		return Lg;
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

	double likelihood(double[][] F, double[][] W) {
		return (1 - alpha) * graphLikelihood(F) + alpha * attributeLikelihood(W, F) - lambda * MathFunctions.L1Norm(W);
	}

	double lineSearchW(HashMap<Integer, ArrayList<Double>> updates) {

		double fx = (alpha) * attributeLikelihood(W, F) - lambda * MathFunctions.L1Norm(W);

		double t = 0.01;
		double beta = 0.3;

		boolean toContinue = true;
		int iter = 0;
		do {
			double[][] W_copy = new double[numAttr][C];
			// System.out.println("t = " + t);
			for (int k = 0; k < numAttr; k++) {
				for (int c = 0; c < C; c++) {
					W_copy[k][c] = W[k][c] + (t * updates.get(k).get(c));
				}
			}

			double fx_new = (alpha) * attributeLikelihood(W_copy, F) - lambda * MathFunctions.L1Norm(W_copy);
			if ((fx_new - fx) < t * MathFunctions.L2NormSq(updates) * 0.01) {
				t = beta * t;
			} else
				toContinue = false;
			iter++;
			if (iter == 10)
				return t;

		} while (toContinue);
		return t;
	}

	double lineSearch(HashMap<Integer, ArrayList<Double>> updates) {
		double beta = 0.6;
		double fx = (1 - alpha) * graphLikelihood(F) + (alpha) * attributeLikelihood(W, F);
		double t = 0.01;
		int iter = 0;
		boolean toContinue = true;
		do {
			double[][] F_copy = new double[V][C];
			boolean nonNegative = true;
			// System.out.println("t = " + t);
			for (int u = 0; u < G.size(); u++) {
				for (int c = 0; c < C; c++) {
					F_copy[u][c] = Math.max(0, F[u][c] + (t * updates.get(u).get(c)));
				}
			}

			if (nonNegative == false) {
				t = beta * t;
				continue;
			}

			double fx_new = (1 - alpha) * graphLikelihood(F_copy) + (alpha) * attributeLikelihood(W, F_copy);
			if ((fx_new - fx) < t * MathFunctions.L2NormSq(updates) * 0.01) {
				t = beta * t;
			} else
				toContinue = false;

			iter++;
			if (iter == 10)
				return t;

		} while (toContinue);
		return t;
	}

	void updateFDriver() throws InterruptedException {
		double[] sumFvc = new double[C];
		HashMap<Integer, ArrayList<Double>> updates = new HashMap<>();

		for (int v = 0; v < G.size(); v++) {
			double[] Fv = F[v];
			for (int c = 0; c < C; c++) {
				sumFvc[c] += Fv[c];
			}
		}
		ArrayList<ArrayList<Integer>> nodeSplits = new ArrayList<>();
		for (int i = 0; i < THREADS; i++)
			nodeSplits.add(new ArrayList<>());

		int i = 0;
		for (int node = 0; node < G.size(); node++) {
			nodeSplits.get(i % THREADS).add(node);
			i++;
		}

		UpdateFThread[] threads = new UpdateFThread[THREADS];

		for (int j = 0; j < threads.length; j++) {
			threads[j] = new UpdateFThread(sumFvc, nodeSplits.get(j));
		}

		for (UpdateFThread thread : threads) {
			thread.start();
		}
		for (UpdateFThread thread : threads) {
			thread.join();
		}

		for (UpdateFThread thread : threads) {
			for (int u : thread.updates.keySet()) {
				updates.put(u, new ArrayList<Double>());
				for (int c = 0; c < C; c++) {
					updates.get(u).add(thread.updates.get(u).get(c));
					// F.get(u).set(c, Math.max(0.0, F.get(u).get(c) + eta *
					// thread.updates.get(u).get(c)));
				}
			}
		}

		// eta = lineSearch(updates);
		System.out.println("Final eta " + eta);
		
		for (int u = 0; u < G.size(); u++) {
			for (int c = 0; c < C; c++) {
				F[u][c] = Math.max(0.0, F[u][c] + eta * updates.get(u).get(c));
			}
		}
	}

	class UpdateFThread extends Thread {

		ArrayList<Integer> nodes;
		double[] sumFvc;
		HashMap<Integer, ArrayList<Double>> updates;

		UpdateFThread(double[] sumFvc, ArrayList<Integer> nodes) {
			this.sumFvc = sumFvc;
			this.nodes = nodes;
		}

		public void run() {
			updates = updateF(sumFvc, nodes);
		}
	}

	class UpdateWThread extends Thread {

		ArrayList<Integer> attributes;
		HashMap<Integer, ArrayList<Double>> updates;

		UpdateWThread(ArrayList<Integer> attributes) {
			this.attributes = attributes;
		}

		public void run() {
			updates = updateW(attributes);
		}
	}

	HashMap<Integer, ArrayList<Double>> updateF(double[] sumFvc, ArrayList<Integer> nodes) {
		// HashMap<Integer, ArrayList<Double>> newF = new HashMap<>();
		HashMap<Integer, ArrayList<Double>> updates = new HashMap<>();
		/*
		 * double[] sumFvc = new double[C];
		 * 
		 * for (int v : G.keySet()) { ArrayList<Double> Fv = F.get(v); for (int
		 * c = 0; c < C; c++) { sumFvc[c] += Fv.get(c); } }
		 */

		for (int u : nodes) {
			// newF.put(u, new ArrayList<Double>());
			updates.put(u, new ArrayList<Double>());
			double[] derivative_G = new double[C];
			double[] derivative_X = new double[C];
			double[] diff = new double[C];
			for (int v : G.get(u)) {
				double dot = MathFunctions.dotProduct(F[u], F[v]);
				double exp = Math.exp(-1 * dot);
				double fraction = exp / (1 - exp);

				for (int c = 0; c < C; c++) {
					derivative_G[c] += F[v][c] * fraction;
					diff[c] -= F[v][c];
					// derivative_G[c] -= F.get(v).get(c);
				}

			}
			for (int c = 0; c < C; c++) {
				diff[c] = diff[c] + sumFvc[c] - F[u][c];
				derivative_G[c] -= diff[c];
			}

			// derivative of Fu done

			for (int k = 0; k < numAttr; k++) {
				int z = X[u][k] ? 1 : 0;
				double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(W[k], F[u])));
				double difference = (z - Quk);
				for (int c = 0; c < C; c++) {
					derivative_X[c] += difference * W[k][c];
				}
			}

			for (int c = 0; c < C; c++) {
				double update = (1 - alpha) * derivative_G[c] + (alpha) * derivative_X[c];
				updates.get(u).add(update);
			}
		}

		// F = new HashMap<>(newF);
		// System.out.println(F);
		return updates;
	}

	// TODO: Function not in use.
	void updateW() {
		HashMap<Integer, ArrayList<Double>> newW = new HashMap<>();

		for (int k = 0; k < numAttr; k++) {
			newW.put(k, new ArrayList<Double>());
			for (int c = 0; c < C; c++) {
				for (int u = 0; u < G.size(); u++) {
					double z = X[u][k] ? 1 : 0;
					double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(W[k], F[u])));
					double diff = (z - Quk);
					double update = 0.0;
					update += diff * F[u][c];
					update *= alpha;
					if (W[k][c] > 0)
						update -= lambda;
					else
						update += lambda;
					newW.get(k).add(update * eta);

				}
			}
		}
		// W = new HashMap<>(newW);
	}

	HashMap<Integer, ArrayList<Double>> updateW(ArrayList<Integer> attributes) {
		HashMap<Integer, ArrayList<Double>> updates = new HashMap<>();
		for (int k : attributes) {
			updates.put(k, new ArrayList<Double>());
			for (int c = 0; c < C; c++) {
				for (int u = 0; u < G.size(); u++) {
					double z = X[u][k] ? 1 : 0;
					double Quk = 1 / (1 + Math.exp(-1 * MathFunctions.dotProduct(W[k], F[u])));
					double diff = (z - Quk);
					double update = 0.0;
					update += diff * F[u][c];
					update *= alpha;
					if (W[k][c] > 0)
						update -= lambda;
					else
						update += lambda;
					updates.get(k).add(update * eta);
				}
			}
		}
		return updates;
	}

	void updateWDriver() throws InterruptedException {
		ArrayList<ArrayList<Integer>> attributeSplits = new ArrayList<>();
		for (int i = 0; i < THREADS; i++)
			attributeSplits.add(new ArrayList<>());

		for (int x = 0; x < numAttr; x++) {
			attributeSplits.get(x % THREADS).add(x);
		}
		UpdateWThread[] threads = new UpdateWThread[THREADS];

		for (int j = 0; j < threads.length; j++) {
			threads[j] = new UpdateWThread(attributeSplits.get(j));
		}

		for (UpdateWThread thread : threads) {
			thread.start();
		}
		for (UpdateWThread thread : threads) {
			thread.join();
		}

		HashMap<Integer, ArrayList<Double>> updates = new HashMap<>();

		for (UpdateWThread thread : threads) {
			for (int k : thread.updates.keySet()) {
				updates.put(k, new ArrayList<Double>());
				for (int c = 0; c < C; c++) {
					updates.get(k).add(thread.updates.get(k).get(c));
				}
			}
		}
		// eta = lineSearchW(updates);
		// eta = 0.0001;
		
		for (int k = 0 ; k < numAttr; k ++) {
			for (int c = 0; c < C; c++) {
				W[k][c] = W[k][c] + eta * updates.get(k).get(c);
			}
		}

	}

	void gradientAscent() throws InterruptedException {
		double startTime = System.currentTimeMillis();

		double previousLikelihood = -1;
		for (int iter = 0; iter < maxIter; iter++) {
			double likelihood = likelihood(F, W);
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
			updateWDriver();
			// updateW();

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
		return 0.0;
		// return likelihood(F, W);
	}

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		main m = new main(10);
		m.driver("facebook/0.edges", "facebook/0.feat");

	}

}
