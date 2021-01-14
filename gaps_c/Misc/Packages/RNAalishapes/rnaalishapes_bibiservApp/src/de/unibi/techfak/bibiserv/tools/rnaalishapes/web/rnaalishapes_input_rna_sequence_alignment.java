/*
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 * 
 * Copyright 2010-2012 BiBiServ Curator Team, http://bibiserv.cebitec.uni-bielefeld.de, 
 * All rights reserved.
 * 
 * The contents of this file are subject to the terms of the Common
 * Development and Distribution License("CDDL") (the "License"). You 
 * may not use this file except in compliance with the License. You can 
 * obtain a copy of the License at http://www.sun.com/cddl/cddl.html
 * 
 * See the License for the specific language governing permissions and 
 * limitations under the License.  When distributing the software, include 
 * this License Header Notice in each file.  If applicable, add the following 
 * below the License Header, with the fields enclosed by brackets [] replaced
 *  by your own identifying information:
 * 
 * "Portions Copyrighted 2010-2012 BiBiServ Curator Team, http://bibiserv.cebitec.uni-bielefeld.de"
 * 
 * Contributor(s):  Daniel Hagemeier - dhagemei[aet]cebitec.uni-bielefeld.de
 *                  Armin Toepfer - atoepfer[aet]cebitec.uni-bielefeld.de
 *                  Jan Krueger - jkrueger[aet]cebitec.uni-bielefeld.de
 * 
 */
package de.unibi.techfak.bibiserv.tools.rnaalishapes.web;


import de.unibi.cebitec.bibiserv.utils.ValidationConnection;
import de.unibi.techfak.bibiserv.web.beans.session.AbstractCloudInputBean;
import de.unibi.cebitec.bibiserv.utils.UniversalRepresentationFinder;
import de.unibi.techfak.bibiserv.util.ontoaccess.bibiontotypes.OntoRepresentation;
import java.util.List;
/**
 * This is a autogenerated bean input template class for function <i>rnaalishapes_input_rna_sequence_alignment</i>
 *
 * @author Daniel Hagemeier - dhagemei[aet]cebitec.uni-bielefeld.de
 *         Armin Toepfer - atoepfer[aet]cebitec.uni-bielefeld.de
 *         Jan Krueger - jkrueger[aet]cebitec.uni-bielefeld.de
 *         Thomas Gatter - tgatter[aet]cebitec.uni-bielefeld.de
 */
public class rnaalishapes_input_rna_sequence_alignment  extends AbstractCloudInputBean {

    
    /** Return id current bean  */
    public String getId(){
        return "rnaalishapes_input_rna_sequence_alignment";
    }
    
    
    /* ########################
     * # validation function  #
     * ######################## */
    public boolean validate(String input, String [] args, OntoRepresentation target){
        List<OntoRepresentation> found = UniversalRepresentationFinder.getOntoRepresentation(input,
                target, this, args);
       if(found.isEmpty()){
           return false;
       }
       
       return true;
    }
    
    public boolean validate(ValidationConnection input, String [] args, OntoRepresentation target){
        List<OntoRepresentation> found = UniversalRepresentationFinder.getOntoRepresentation(input,
                target, this, args);
       if(found.isEmpty()){
           return false;
       }
       
       return true;
    }  
    
    
    /**
     * TODO: TG 01/13 Change this if stream support changes.
     * @return true if the tool supports streams as input (stdin or unix named pipe file)
     */
    public boolean supportsStreamedInput() {
        return false;
    }
}
