package org.jlab.clas12.fiducials;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import javax.swing.JFrame;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.base.GeometryFactory;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.geom.base.ConstantProvider;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Trap3D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.*;
import org.jlab.logging.DefaultLogger;
import org.jlab.utils.benchmark.ProgressPrintout;
import org.jlab.utils.options.OptionParser;
        
        
/**
 *
 * @author devita
 */
public class Fiducial {
    
    Map<Integer, Map<Integer,DataGroup>> dgs = new LinkedHashMap<>();
    
    double delta = 3;
    double[] thetas = {8, 15, 23};
    String[] types = {"Data", "MC"};

    DCGeant4Factory dcFactory = null;
    
    Random rand = new Random();
    
    public Fiducial() {
        this.createHistos();
        
        ConstantProvider dcProvider   = GeometryFactory.getConstants(DetectorType.DC, 11, "rga_fall2018");
        dcFactory = new DCGeant4Factory(dcProvider, DCGeant4Factory.MINISTAGGERON, false);

    }
    
    
    
    
    public final void createHistos() {

        GStyle.getAxisAttributesX().setTitleFontSize(24); //24
        GStyle.getAxisAttributesX().setLabelFontSize(18); //18
        GStyle.getAxisAttributesY().setTitleFontSize(24); //24
        GStyle.getAxisAttributesY().setLabelFontSize(18); //18
        GStyle.getAxisAttributesZ().setLabelFontSize(14); //14
        GStyle.setPalette("kDefault");
        GStyle.getAxisAttributesX().setLabelFontName("Arial");
        GStyle.getAxisAttributesY().setLabelFontName("Arial");
        GStyle.getAxisAttributesZ().setLabelFontName("Arial");
        GStyle.getAxisAttributesX().setTitleFontName("Arial");
        GStyle.getAxisAttributesY().setTitleFontName("Arial");
        GStyle.getAxisAttributesZ().setTitleFontName("Arial");
        GStyle.setGraphicsFrameLineWidth(1);
        GStyle.getH1FAttributes().setLineWidth(2);
        GStyle.getH1FAttributes().setOptStat("");

        DataGroup dgEB   = new DataGroup(2, 2);
        DataGroup dgTraj = new DataGroup(3, 2);
        DataGroup dgPhi  = new DataGroup(3, 2);
        DataGroup dgEdge = new DataGroup(3, 2);
        DataGroup dgPhiT = new DataGroup(3, 2);
        DataGroup dgEdgT = new DataGroup(3, 2);
        DataGroup dgCutT = new DataGroup(3, 2);
        for(int type=0; type<2; type++) {
            H2F hiptheta   = new H2F("hiptheta"+type, types[type], 200, 0, 6, 200, 5, 45);
            hiptheta.setTitleX("p (GeV)");
            hiptheta.setTitleY("#theta (deg)");
            H2F hiphitheta = new H2F("hiphitheta"+type, types[type], 200, -180, 180, 200, 5, 45);
            hiphitheta.setTitleX("#phi (deg)");
            hiphitheta.setTitleY("#theta (deg)");
            dgEB.addDataSet(hiptheta,   0+type*2);
            dgEB.addDataSet(hiphitheta, 1+type*2);
            for(int i=0; i<3; i++) {
                int region = i+1;
                double size = 90+i*60;
                H2F hi1 = new H2F("hi"+type+region, "Region "+region, 200, -size, size, 200, -size, size);  
                hi1.setTitleX("y (cm)"); 
                hi1.setTitleY("x (cm)"); 
                H2F hi2 = new H2F("hi"+type+region, "Region"+region, 200, -30, 30, 200, 0, 50);  
                hi2.setTitleX("#phi_t_r_a_j (deg)"); 
                hi2.setTitleY("#theta_t_r_a_j (deg)"); 
                H2F hi3 = new H2F("hi"+type+region, "Region"+region, 200, -size, size, 200, -size, size);  
                hi3.setTitleX("y (cm)"); 
                hi3.setTitleY("x (cm)"); 
                for(int k=0; k<4; k++) {
                    H1F hi4 = new H1F("hi"+type+region+k, "Region"+region, 100, -30, 30);  
                    hi4.setTitleX("#phi_t_r_a_j (deg)"); 
                    hi4.setTitleY("Counts"); 
                    hi4.setLineColor(k+1);
                    dgPhiT.addDataSet(hi4, i+type*3);
                }
                for(int k=0; k<2; k++) {
                    H1F hi5 = new H1F("hi"+type+region+k, "Region"+region, 40, -5, 15);  
                    hi5.setTitleX("edge (cm)"); 
                    hi5.setTitleY("Counts"); 
                    hi5.setLineColor(k+1);
                    dgEdgT.addDataSet(hi5, i+type*3);
                }
                for(int k=0; k<2; k++) {
                    H1F hi6 = new H1F("hi"+type+region+k, "R1 - #theta_t_r_a_j = "+thetas[i]+" (deg)", 40, -5, 15);  
                    hi6.setTitleX("edge (cm)"); 
                    hi6.setTitleY("Counts"); 
                    hi6.setLineColor(k+1);
                    dgCutT.addDataSet(hi6, i+type*3);
                }
                dgTraj.addDataSet(hi1, i+type*3);
                dgPhi.addDataSet(hi2,  i+type*3);
                dgEdge.addDataSet(hi3, i+type*3);
            }
        }
        if(!dgs.containsKey(DetectorType.HTCC.getDetectorId()))
            dgs.put(DetectorType.HTCC.getDetectorId(), new LinkedHashMap());
        dgs.get(DetectorType.HTCC.getDetectorId()).put(1, dgEB);
        if(!dgs.containsKey(DetectorType.DC.getDetectorId()))
            dgs.put(DetectorType.DC.getDetectorId(), new LinkedHashMap());
        dgs.get(DetectorType.DC.getDetectorId()).put(1, dgTraj);
        dgs.get(DetectorType.DC.getDetectorId()).put(2, dgEdge);
        dgs.get(DetectorType.DC.getDetectorId()).put(3, dgPhi);
        dgs.get(DetectorType.DC.getDetectorId()).put(4, dgPhiT);
        dgs.get(DetectorType.DC.getDetectorId()).put(5, dgEdgT);
        dgs.get(DetectorType.DC.getDetectorId()).put(6, dgCutT);
    }
    
    
    public void processEvent(DataEvent event, int type) {
        DataBank recPart = null;
        DataBank recTraj = null;
        DataBank recTrac = null;
        DataBank mcPart = null;
        DataBank mcTrue = null;

        if (!event.hasBank("RUN::config")) {
            return;
        }
        int run = event.getBank("RUN::config").getInt("run", 0);

        if (event.hasBank("MC::True")) {
            mcTrue = event.getBank("MC::True");
        }
        if (event.hasBank("MC::Particle")) {
            mcPart = event.getBank("MC::Particle");
        }
        if (event.hasBank("REC::Particle")) {
            recPart = event.getBank("REC::Particle");
        }
        if (event.hasBank("REC::Traj")) {
            recTraj = event.getBank("REC::Traj");
        }
        if (event.hasBank("REC::Track")) {
            recTrac = event.getBank("REC::Track");
        }
        
        if(mcPart!=null) {
            
        }
        if (recPart != null && recTraj!=null && recTrac!=null) {
//            if(recPart.getInt("pid", 0)!=11 || recPart.getShort("status", 0)>-2000) return;
            for(int i=0; i<recTraj.rows(); i++) {
                int pindex   = recTraj.getShort("pindex", i);
                int pid      = recPart.getInt("pid", pindex);
                int status   = recPart.getShort("status", pindex);
                int clas     = Math.abs(status)/1000;
                if(pid!=211 || clas!=2) continue;
                int detector = recTraj.getInt("detector", i);
                int layer    = recTraj.getInt("layer", i);
                double x     = recTraj.getFloat("x", i);
                double y     = recTraj.getFloat("y", i);
                double z     = recTraj.getFloat("z", i);
                double edge  = recTraj.getFloat("edge", i);
                int sector = 0;
                for(int j=0; j<recTrac.rows(); j++) {
                    if(recTrac.getShort("pindex", j)==pindex)
                        sector = recTrac.getByte("sector", j);
                }
                if(sector!=2) continue;
                Particle part = new Particle(pid,
                                             recPart.getFloat("px", pindex),
                                             recPart.getFloat("py", pindex),
                                             recPart.getFloat("pz", pindex),
                                             recPart.getFloat("vx", pindex),
                                             recPart.getFloat("vy", pindex),
                                             recPart.getFloat("vz", pindex));
                
                Point3D trajTilted = new Point3D(x, y, z);
                Point3D trajLocal  = new Point3D(x, y, z);
                trajTilted.rotateZ(-Math.PI/3*(sector-1));
                trajLocal.rotateZ(Math.PI/2-Math.PI/3*(sector-1));
                if(layer==3 && detector==DetectorType.FTOF.getDetectorId())
                    trajLocal.rotateX(Math.toRadians(58.11));
                else
                    trajLocal.rotateX(Math.toRadians(25));
                
                if(detector==DetectorType.DC.getDetectorId() && (layer==6 || layer==18 || layer==36)) {
                    int region = (layer-1)/12+1;
                    double theta = Math.toDegrees(trajTilted.toVector3D().theta());
                    double phi   = Math.toDegrees(trajTilted.toVector3D().phi());
                    //System.out.println(region);
                    dgs.get(detector).get(1).getH2F("hi"+type+region).fill(trajLocal.x(), trajLocal.y());
                    if(edge>delta) dgs.get(detector).get(2).getH2F("hi"+type+region).fill(trajLocal.x(), trajLocal.y());
                    dgs.get(detector).get(3).getH2F("hi"+type+region).fill(phi, theta);
                    if(Math.abs(theta-15)<1){
                        for(int k=0; k<4; k++) {
                            if(k==0 || edge>-delta+delta*k) 
                                dgs.get(detector).get(4).getH1F("hi"+type+region+k).fill(phi);
                        }
                        if(edge>0) {
                            int side = (int) (1-Math.signum(trajLocal.x()))/2;
                            dgs.get(detector).get(5).getH1F("hi"+type+region+side).fill(edge);
                        }
                    }
                    for(int r=0; r<3; r++) {
                        if(edge>0 && Math.abs(theta-(8+7*r))<1) {
                            int side = (int) (1-Math.signum(trajLocal.x()))/2;
                            dgs.get(detector).get(6).getH1F("hi"+type+(r+1)+side).fill(edge);
                        }
                    }
                }
                else if(detector==DetectorType.HTCC.getDetectorId()) {
                    dgs.get(detector).get(1).getH2F("hiptheta"+type).fill(part.p(), Math.toDegrees(part.theta()));
                    dgs.get(detector).get(1).getH2F("hiphitheta"+type).fill(Math.toDegrees(part.phi()), Math.toDegrees(part.theta()));
                }
            }		
        }

    }
    
    
    public void analyzeHisto(int nevents) {
        
        for(int type=0; type<2; type++) {
            for(int region=1; region<4; region++) {
                for(int side=0; side<2; side++) {
                    H1F hi = dgs.get(DetectorType.DC.getDetectorId()).get(5).getH1F("hi"+type+region+side);
    //                hi.divide((hi.getBinContent(hi.getDataSize(0)-1)+hi.getBinContent(hi.getDataSize(0)-2)+hi.getBinContent(hi.getDataSize(0)-3))/3);
                }
            }
        }
//        nevents *= current/50;
//
//        double norm = nevents*6/100;
//        for(int i=0; i<titles.length; i++) {
//            String title = titles[i];
//
//            H2F occ = dgs.get("DC").getH2F("hiOcc2D"+title);
//            for(int loop = 0; loop < occ.getDataBufferSize(); loop++){
//                double r = occ.getDataBufferBin(loop);
//                if(r>0) occ.setDataBufferBin(loop,r/norm);
//            }
//
//            dgs.get("DC").getH1F("hiOcc1D"+title).divide((double) nevents*6*12*112/100);
//
//            for(int iregion=0; iregion<3; iregion++) {
//                int region = iregion +1;
//                dgs.get("DC").getH1F("hiTDC" + title + region).divide((double) nevents);
//            }
//        }
//        H2F occ = dgs.get("URWELL").getH2F("hiOcc2D");
//        for(int loop = 0; loop < occ.getDataBufferSize(); loop++){
//            double r = occ.getDataBufferBin(loop);
//            if(r>0) occ.setDataBufferBin(loop,r/norm);
//        }
    }
    
        
    public EmbeddedCanvasTabbed getCanvas() {
        EmbeddedCanvasTabbed canvas = null;
        for(int key : dgs.keySet()) {
            for(int tab : dgs.get(key).keySet()) {
                if(canvas == null)
                    canvas = new EmbeddedCanvasTabbed(DetectorType.getType(key).getName()+tab);
                else
                    canvas.addCanvas(DetectorType.getType(key).getName()+tab);
                canvas.getCanvas(DetectorType.getType(key).getName()+tab).draw(dgs.get(key).get(tab));
                for(EmbeddedPad pad : canvas.getCanvas(DetectorType.getType(key).getName()+tab).getCanvasPads()) {
                    pad.setTitleFont("Arial");
                    pad.setTitleFontSize(18);
                    if(pad.getDatasetPlotters().size()>0 && pad.getDatasetPlotters().get(0).getDataSet() instanceof H2F) 
                        pad.getAxisZ().setLog(true);
                }
                if(key==DetectorType.DC.getDetectorId() && tab<3) this.drawDC(canvas.getCanvas(DetectorType.getType(key).getName()+tab));
            }
        }
        return canvas;
    }
    
    private void drawDC(EmbeddedCanvas canvas) {
        
        for(int layer : new int[]{6, 18, 36}) {
            int superlayer = layer/6;
            for(int iedge=0; iedge<2; iedge++) {
                canvas.cd((superlayer-1)/2+iedge*3);
                Trap3D contour = dcFactory.getTrajectorySurface(0, superlayer-1, (layer-1)%6);
                contour.rotateZ(Math.PI/2);
                contour.rotateX(Math.toRadians(25));
                DataLine[] lines = new DataLine[4];
                for(int i=0; i<lines.length; i++) {
                    int i1 = i;
                    int i2 = i+1; if(i2==lines.length) i2=0;
                    lines[i] = new DataLine(contour.point(i1).x(),contour.point(i1).y(),contour.point(i2).x(),contour.point(i2).y());
                    lines[i].setLineColor(2);
                    canvas.draw(lines[i]);
                }
            }
        }
    }

    public static void main(String[] args) {
        

        OptionParser parser = new OptionParser("fiducials");
        parser.setRequiresInputList(true);
        parser.addOption("-n","-1", "maximum number of events to process");
//        parser.addOption("-w", "1", "open graphical window (1) or run in batch mode (0)");
        parser.parse(args);
        
        
        int     maxEvents = parser.getOption("-n").intValue();
        boolean window    = true;//parser.getOption("-w").intValue()==1;
        
        if(!window) System.setProperty("java.awt.headless", "true");
        DefaultLogger.debug();
        
        Fiducial analysis = new Fiducial();        

        List<String> inputFiles = parser.getInputList();
        
        ProgressPrintout progress = new ProgressPrintout();

        int counter=-1;

        for(int type=0; type<inputFiles.size(); type++) {
            String input = inputFiles.get(type);
            HipoDataSource  reader = new HipoDataSource();
            reader.open(input);
            while(reader.hasEvent()) {
                
                counter++;
                
                DataEvent ev = reader.getNextEvent();

                analysis.processEvent(ev, type);
                
                progress.updateStatus();
                
                if(maxEvents>0){
                    if(counter>=maxEvents) break;
                }

            }
            progress.showStatus();
            reader.close();
        }   
        analysis.analyzeHisto(counter);
    
        if(window) {
            JFrame frame = new JFrame("Fiducials");
            frame.setSize(1500,1000);
            frame.add(analysis.getCanvas());
            frame.setLocationRelativeTo(null);
            frame.setVisible(true);     
        }
    }
}



