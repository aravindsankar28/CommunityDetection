
public class UpdateWThread extends Thread {
	double[] updates;
	int startIndex;
	int endIndex;
	double[][] F;
	double[] W;

	UpdateWThread(int startIndex, int endIndex, double[] updates, double[][] F, double[] W) {
		this.startIndex = startIndex;
		this.endIndex = endIndex;
		this.updates = updates;
		this.F = F;
		this.W = W;
	}

	public void run() {
		updates = updateW(startIndex, endIndex, updates, F, W);
	}

	// The gradient updates are computed here for every W_c for all the
	// communities assigned to this thread as per the equation described.

	double[] updateW(int startIndex, int endIndex, double[] updates, double[][] F, double[] W) {
		for (int c = startIndex; c <= endIndex; c++) {
			double update = 0.0;
			for (int u = 0; u < Model.G.size(); u++) {
				for (int v = u + 1; v < Model.G.size(); v++) {
					if (Model.G.get(u).containsKey(v)) {
						double psi = MathFunctions.computePsi(u, v, F[u], F[v], W, Model.cosineSim, Model.alpha);
						double exp = Math.exp(-1 * psi);
						double fraction = exp / (1 - exp);
						update += fraction * F[u][c] * F[v][c];
					} else
						update -= F[u][c] * F[v][c];
				}
			}
			update *= Model.alpha;
			updates[c] = update;
		}
		return updates;
	}

}
