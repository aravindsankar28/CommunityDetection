import java.util.ArrayList;

public class CesnaMain {
	static main m;
	int minCommunities;
	int maxCommunities;
	static String graphfilename;
	static String attrfilename;

	CesnaMain() {
		graphfilename = "facebook/0.edges";
		attrfilename = "facebook/0.feat";

		minCommunities = 5;
		maxCommunities = 15;
	}

	int estimateCommunitiesBIC() throws Exception {
		ArrayList<Double> BIC_values = new ArrayList();
		for (int i = minCommunities; i < maxCommunities; i++) {
			m = new main(i);
			double likelihood = m.driver(graphfilename, attrfilename);
			int k = m.numberOfParameters();
			double BIC = -2* likelihood + k* Math.log(m.V);
			System.out.println(BIC);
			BIC_values.add(BIC);
		}
		System.out.println(BIC_values);
		
		return 0;
	}

	public static void main(String[] args) throws Exception {
		CesnaMain object = new CesnaMain();
		object.estimateCommunitiesBIC();
		// m.driver(graphfilename, attrfilename);
	}

}
