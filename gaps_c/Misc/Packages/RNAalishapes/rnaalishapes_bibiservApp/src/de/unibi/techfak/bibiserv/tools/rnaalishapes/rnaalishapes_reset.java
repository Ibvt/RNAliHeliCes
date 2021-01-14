/*
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Copyright 2010 BiBiServ Curator Team, http://bibiserv.cebitec.uni-bielefeld.de,
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
 * "Portions Copyrighted [year] [name of copyright owner]"
 *
 * Contributor(s):
 *
 */

package de.unibi.techfak.bibiserv.tools.rnaalishapes;

import java.io.IOException;
import java.util.Map;
import org.apache.log4j.Logger;
import javax.faces.context.FacesContext;

/**
 * A bean used to reset the session for just this tool.
 * @author Thomas Gatter - tgatter(at)cebitec.uni-bielefeld.de
 */
public class rnaalishapes_reset {

    private static Logger log = Logger.getLogger(rnaalishapes_reset.class);


    public void reset() {
        
       String[] beanIds = {"awsBean", "toolBean_rnaalishapes_download", "rnaalishapes_input_rna_sequence_alignment", "rnaalishapes_input_rna_secondary_structure", "rnaalishapes_function", "rnaalishapes_function_mfe_param", "rnaalishapes_function_mfe_result", "rnaalishapes_function_mfe_resulthandler", "rnaalishapes_function_mfe_execfunction", "rnaalishapes_function_mfe_controller", "rnaalishapes_function_subopt_param", "rnaalishapes_function_subopt_result", "rnaalishapes_function_subopt_resulthandler", "rnaalishapes_function_subopt_execfunction", "rnaalishapes_function_subopt_controller", "rnaalishapes_function_shapes_param", "rnaalishapes_function_shapes_result", "rnaalishapes_function_shapes_resulthandler", "rnaalishapes_function_shapes_execfunction", "rnaalishapes_function_shapes_controller", "rnaalishapes_function_probs_param", "rnaalishapes_function_probs_result", "rnaalishapes_function_probs_resulthandler", "rnaalishapes_function_probs_execfunction", "rnaalishapes_function_probs_controller", "rnaalishapes_function_sample_param", "rnaalishapes_function_sample_result", "rnaalishapes_function_sample_resulthandler", "rnaalishapes_function_sample_execfunction", "rnaalishapes_function_sample_controller", "rnaalishapes_function_eval_param", "rnaalishapes_function_eval_result", "rnaalishapes_function_eval_resulthandler", "rnaalishapes_function_eval_execfunction", "rnaalishapes_function_eval_controller", "rnaalishapes_function_convert_param", "rnaalishapes_function_convert_result", "rnaalishapes_function_convert_resulthandler", "rnaalishapes_function_convert_execfunction", "rnaalishapes_function_convert_controller", "rnaalishapes_function_outside_param", "rnaalishapes_function_outside_result", "rnaalishapes_function_outside_resulthandler", "rnaalishapes_function_outside_execfunction", "rnaalishapes_function_outside_controller"};
        
        Map<String,Object> sessionMap = FacesContext.getCurrentInstance().getExternalContext().getSessionMap();
        
        for(String id: beanIds) {
            if(sessionMap.containsKey(id)) {
                sessionMap.remove(id);
            }
        }
        
        try {
            FacesContext.getCurrentInstance().getExternalContext().redirect("/rnaalishapes");
        } catch (IOException ex) {
            log.warn("Could not redirect on reset: "+ex);
        }
    }

}
