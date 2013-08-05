import java.io.*;
import java.util.*;

public class MatlabCAPRIParser {
	public static void main(String[] args) throws IOException {
		String prot = args[0];
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("results\\capri\\result" + prot + "\\result" + prot + "matlab.txt")));
		BufferedReader br1 = new BufferedReader(new FileReader("results\\capri\\result" + prot.toUpperCase() + "\\scorerms" + prot + ".txt"));
		double[] rms = new double[10];
		double[] score = new double[10];
		for (int i = 0; i < 10; i++) {
			StringTokenizer st = new StringTokenizer(br1.readLine());
			score[i] = Double.parseDouble(st.nextToken());
			rms[i] = Double.parseDouble(st.nextToken());
		}
		out.print("x1=[");
		for (int i = 0; i < 10; i++) {
			out.print(rms[i]);
			if (i != 9) out.print(",");
		}
		out.print("];");
		out.println();
		out.print("y1=[");
		for (int i = 0; i < 10; i++) {
			double test = score[i];
			if (!Double.isInfinite(test) && !Double.isNaN(test)) out.print(test);
			else out.print(10000);
			if (i != 9) out.print(",");
		}
		out.print("];");
		out.println();
		BufferedReader br2 = new BufferedReader(new FileReader("results\\capri\\result" + prot.toUpperCase() + "\\scorermsiter" + prot + ".txt"));
		double[] rms2 = new double[10];
		for (int i = 0; i < 10; i++) {
			StringTokenizer st = new StringTokenizer(br2.readLine());
			st.nextToken();
			rms2[i] = Double.parseDouble(st.nextToken());
		}
		out.print("x2=[");
		for (int i = 0; i < 10; i++) {
			out.print(rms2[i]);
			if (i != 9) out.print(",");
		}
		out.print("];");
		out.println();
		out.print("y2=[");
		for (int i = 0; i < 10; i++) {
			double test = score[i];
			if (!Double.isInfinite(test) && !Double.isNaN(test)) out.print(test);
			else out.print(10000);
			if (i != 9) out.print(",");
		}
		out.print("];");
		out.println();
		out.println("figure1 = figure; \naxes1 = axes('Parent',figure1,'LineWidth',5,'FontWeight','bold','FontSize',36,'FontName','TimesNewRoman','CLim',[0 1]); \nhold(axes1,'all'); \nscatter(x1,y1,100,'black','MarkerFaceColor','flat','MarkerEdgeColor','none','Parent',axes1,'DisplayName','Without iteration'); \nscatter(x2,y2,100,'red','MarkerFaceColor','flat','MarkerEdgeColor','none','Parent',axes1,'DisplayName','With iteration'); \nxlabel({'LRMS (A)'},'FontWeight','bold','FontSize',72,'FontName','Times New Roman'); \nylabel({'Score'},'FontWeight','bold','FontSize',72,'FontName','Times New Roman'); \ntitle({'CAPRI Target " + prot + "'},'FontWeight','bold','FontSize',72,'FontName','Times New Roman');");
		out.flush();
		out.close();
	}
}
