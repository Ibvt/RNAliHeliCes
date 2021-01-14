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
 * Download Item class stores information about a download item like
 * name, filename, description, size (in kb) and (mime-)type. name and
 * filename.
 *
 *
 * @author Jan Krueger - jkrueger(at)cebitec.uni-bielefeld.de
 */
public class DownloadItem {

    /* id of download item  - optional*/
    private String id = null;
    /* name of download item */
    private String name;
    /* filename of download item */
    private String filename;
    /* description of download item - optional*/
    private String description;
    /* download size (in kb) of download item */
    private int size;
    /* (mime) type of download item */
    private String type;


    public DownloadItem(String name, String filename, String description){
        this(name,filename,description,-1,"unknown");
    }

    public DownloadItem(String name, String filename, String description, int size, String type){
        this.name = name;
        this.filename = filename;
        this.description = description;
        this.size = size;
        this.type = type;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getFilename() {
        return filename;
    }

    public void setFilename(String filename) {
        this.filename = filename;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getSize() {
        return size;
    }

    public void setSize(int size) {
        this.size = size;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }


}
