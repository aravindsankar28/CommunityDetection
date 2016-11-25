import java.util.ArrayList;
import java.util.HashMap;

public class MathFunctions {

	static double dotProduct(ArrayList<Double> x, ArrayList<Double> y) {
		double sum = 0.0;
		if (x.size() != y.size())
			return -1;

		for (int i = 0; i < x.size(); i++) {
			sum += x.get(i) * y.get(i);
		}
		return sum;
	}

	static double dotProduct(double[] x, double[] y) {
		double sum = 0.0;
		if (x.length != y.length)
			{
				System.out.println("error");
				return -1;
			}

		for (int i = 0; i < x.length; i++) {
			sum += x[i] * y[i];
		}
		return sum;
	}

	static double L1Norm(HashMap<Integer, ArrayList<Double>> M) {
		double s = 0.0;
		for (int x : M.keySet()) {
			for (double y : M.get(x))
				s += Math.abs(y);
		}
		return s;
	}

	static double L1Norm(double[][] M) {
		double s = 0.0;
		for (int x = 0; x < M.length; x++) {
			for (int y = 0; y < M[x].length; y++)
				s += Math.abs(M[x][y]);
		}
		return s;
	}

	static double L2NormSq(HashMap<Integer, ArrayList<Double>> M) {
		double s = 0.0;
		for (int x : M.keySet()) {
			for (double y : M.get(x))
				s += y * y;
		}
		return s;
	}

	static double L2NormSq(double[][] M) {
		double s = 0.0;
		for (int x = 0; x < M.length; x++) {
			for (int y = 0; y < M[x].length; y++) {
				s += M[x][y] * M[x][y];
			}
		}
		return s;
	}
}
