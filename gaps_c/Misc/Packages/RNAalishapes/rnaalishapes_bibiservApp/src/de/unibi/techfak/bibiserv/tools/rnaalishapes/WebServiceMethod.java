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


/**
 * WebServiceMethod Backing Bean class.
 *
 * Attention! Methods isn't finished, yet! Need access to WS Endpoint implementation
 * object to obtain supported ws methods.
 * 
 *
 *
 * @author Jan Krueger - jkrueger[aet]cebitec.uni-bielefeld.de
 */
public class WebServiceMethod {

    public String[]  methodList(String ftc_id){
       return new String [] {
            "Functionality not yet implemented!",
            "Need access to WSEndpoint bean/object to",
            "determine all public accessible methods."
        };
    }

    public String wsdl(String fct_id){
        return  "/services/ws_"+fct_id+"?wsdl";
    }


    /* #########################################
     * #     DI method setToolDescription      #
     * ######################################### */
    private Tooldescription tooldescription;

    public void setTooldescription(Tooldescription tooldescription) {
        this.tooldescription = tooldescription;
    }

}
