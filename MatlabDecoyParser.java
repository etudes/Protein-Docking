import java.io.*;
import java.util.*;

public class MatlabDecoyParser {
	public static void main(String[] args) throws IOException {
		String prot = args[0];
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("results\\decoys\\" + prot + "matlab.txt")));
		BufferedReader br1 = new BufferedReader(new FileReader("results\\decoys\\" + prot.toUpperCase() + "decoys.txt"));
		BufferedReader br2 = new BufferedReader(new FileReader("results\\decoys\\aa" + prot.toUpperCase() + ".fasc"));
		double[] rms = new double[1000];
		while (true) {
			String s = br2.readLine();
			if (s == null) break;
			StringTokenizer st = new StringTokenizer(s);
			String test = st.nextToken();
			if (!test.startsWith("aa")) continue;
			int num = Integer.parseInt(test.substring(11, 15));
			st.nextToken();
			rms[num-1] = Double.parseDouble(st.nextToken());
		}
		out.print("x1=[");
		for (int i = 0; i < 1000; i++) {
			out.print(rms[i]);
			if (i != 999) out.print(",");
		}
		out.print("];");
		out.println();
		out.print("y1=[");
		for (int i = 0; i < 1000; i++) {
			out.print(Double.parseDouble(br1.readLine()));
			if (i != 999) out.print(",");
		}
		out.print("];");
		out.println();
		BufferedReader br3 = new BufferedReader(new FileReader("decoySummary.txt"));
		boolean failed = false;
		ArrayList<Integer> runKB = new ArrayList<Integer>();
		while (true) {
			String s = br3.readLine();
			if (s == null) break;
			StringTokenizer st = new StringTokenizer(s);
			st.nextToken();
			String test = st.nextToken();
			if (!test.substring(0, 4).toLowerCase().equals(prot)) continue;
			while (st.hasMoreTokens()) {
				test = st.nextToken();
				if (test.equals("SUCCESS")) { failed = true; break; }      
				if (!test.equals("FAILED")) runKB.add(Integer.parseInt(test));
				if (test.equals("FAILED")) break;
			}
		}
		BufferedReader br4 = new BufferedReader(new FileReader("results\\decoys\\" + prot.toUpperCase() + "decoyskb.txt"));
		out.print("y2=[");
		boolean first = true;
		while (true) {
			String s = br4.readLine();
			if (s == null) break;
			if (first) first = false;
			else out.print(",");
			out.print(Double.parseDouble(s));
		}
		out.print("];");
		out.println();
		out.print("x2=[");
		for (int i = 0; i < runKB.size(); i++) {
			out.print(rms[runKB.get(i)]);
			if (i != runKB.size() - 1) out.print(",");
		}
		out.print("];");
		out.println();
		BufferedReader br5 = new BufferedReader(new FileReader("results\\natives\\" + prot + "native.txt"));
		out.println("x3=[0]");
		out.println("y3=[" + Double.parseDouble(br5.readLine()) + "]");
		out.println("figure1 = figure; \naxes1 = axes('Parent',figure1,'CLim',[0 1]); \nhold(axes1,'all'); \nscatter(x1,y1,10,'black','MarkerFaceColor','flat','MarkerEdgeColor','none','Parent',axes1,'DisplayName','Without iteration'); \nscatter(x2,y2,10,'red','MarkerFaceColor','flat','MarkerEdgeColor','none','Parent',axes1,'DisplayName','With iteration'); \nscatter(x3,y3,10,'blue','MarkerFaceColor','flat','MarkerEdgeColor','none','Parent',axes1,'DisplayName','Native configuration'); \nxlabel({'LRMS (A)'},'FontSize',36,'FontName','Times New Roman','Color',[0 0.498039215686275 0]); \nylabel({'Score'},'FontSize',36,'FontName','Times New Roman','Color',[0 0.498039215686275 0]); \ntitle({'Complex " + prot + " decoys'},'FontSize',36,'FontName','Times New Roman','Color',[0 0.498039215686275 0]); \nlegend(axes1,'show');");
		out.flush();
		out.close();
	}
}
