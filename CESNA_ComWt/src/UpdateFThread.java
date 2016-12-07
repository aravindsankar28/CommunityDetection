public class UpdateFThread extends Thread {
	double[] sumFvc;
	int startIndex;
	int endIndex;
	double[][] updates;
	double[][] F;
	double[] W;

	UpdateFThread(double[] sumFvc, int startIndex, int endIndex, double[][] updates, double[][] F, double[] W) {
		this.sumFvc = sumFvc;
		this.startIndex = startIndex;
		this.endIndex = endIndex;
		this.updates = updates;
		this.F = F;
		this.W = W;
	}

	public void run() {
		updates = updateF(sumFvc, startIndex, endIndex, updates, F, W);
	}

	// The gradient updates are computed here for every F_uc for all the users
	// assigned to this thread as per the equation described.
	double[][] updateF(double[] sumFvc, int startIndex, int endIndex, double[][] updates, double[][] F, double[] W) {
		int C = sumFvc.length;
		double alpha = Model.alpha;
		for (int u = startIndex; u <= endIndex; u++) {
			double[] derivative_G = new double[C];
			double[] diff = new double[C];
			for (int v : Model.G.get(u).keySet()) {

				double psi = MathFunctions.computePsi(u, v, F[u], F[v], W, Model.cosineSim, alpha);

				double exp = Math.exp(-1 * psi);
				double fraction = exp / (1 - exp);

				for (int c = 0; c < C; c++) {
					derivative_G[c] += F[v][c] * fraction * alpha * W[c];
					diff[c] += F[v][c];
				}
			}
			for (int c = 0; c < C; c++) {
				derivative_G[c] -= (sumFvc[c] - diff[c] - F[u][c]) * alpha * W[c];
			}

			for (int c = 0; c < C; c++) {
				double update = derivative_G[c];
				updates[u][c] = update;
			}
		}
		return updates;
	}
}