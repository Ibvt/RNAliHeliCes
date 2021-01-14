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

import de.unibi.techfak.bibiserv.cms.Tfile;
import de.unibi.techfak.bibiserv.web.beans.session.MessagesInterface;
import java.util.List;
import org.apache.log4j.Logger;

/**
 * Download Item class stores information about a download item like
 * name, filename, description, size (in kb) and (mime-)type. name and
 * filename.
 *
 *
 * @author Jan Krueger - jkrueger(at)cebitec.uni-bielefeld.de
 */
public class DownloadList {

    private static Logger log = Logger.getLogger(DownloadList.class);


    /**
     * Return an array of DownloadItems. This list can't be hardcoded because
     * of i18n messages ...
     *
     * @return a array of DownloadItems
     */
    public DownloadItem[] getDownloadItems() {
        return generateDownloadItems();
    }


    /* #########################################
     * #     DI method setMessages             #
     * ######################################### */
    private MessagesInterface messages;

    public void setMessages(MessagesInterface messages) {
        this.messages = messages;
    }

    /* #########################################
     * #     DI method setToolDescription      #
     * ######################################### */
    private Tooldescription tooldescription;

    public void setTooldescription(Tooldescription tooldescription) {
        this.tooldescription = tooldescription;
    }


   /* #########################################
    * #        private helper method(s)       #
    * #########################################
    */

    private DownloadItem[] generateDownloadItems(){

        List<Tfile> list_of_downloadable = tooldescription.getToolDescription().getDownloadable();

        DownloadItem [] downloaditems = new DownloadItem [list_of_downloadable.size()];

        for (int i = 0; i < list_of_downloadable.size(); ++i){
            Tfile downloadable = list_of_downloadable.get(i);
            String id = downloadable.getId();
            downloaditems[i] = new DownloadItem(messages.property(id+"_name"),
                                                "/applications/"+tooldescription.getToolDescription().getId()+"/resources/downloads/"+downloadable.getFilename(),
                                                messages.property(id+"_shortDescription"));
        }
        return downloaditems;
    }



}
