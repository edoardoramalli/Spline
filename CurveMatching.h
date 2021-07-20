class CurveMatching {

public:

    /* Runs the Curve Matching program */
    void run();

////////////////////////////////////////////////////////////////////////////////

private:

    /* List of the paths to the folders containing the input files */
    vector<string> folderPathList;

    /* List of the folders containing the input files */
    vector<string> folderList;

    /* List of the names of the files contained in each folder of input data,
    minus the _exp.txt and _mod.txt appendixes. For files 'Example_exp.txt' and
    'Example_mod.txt', nameLists would contain the string 'Example'.
    nameLists[i] refers to folderList[i] */
    vector<vector<string>> nameLists;

    /* Indexes calculated for each pair of input files (_exp and _mod).
    indexes[i] refers to folderList[i] */
    vector<vector<Indexes>> indexes;

    /* Specifies whether the shift index is to be calculated for a specific
    folder. Ordered the same way as folderList */
    vector<bool> calculateShift;

    /* Date when the Curve Matching program starts running */
    string date;

    /* Time when the Curve Matching program starts running */
    string time;

    ////////////////////////////////////////////////////////////////////////////

    
    /* Obtains the list of the non-invisible folders at a specific path and the
    list of their paths */
    void listFoldersListFolderPaths(const string& path);

    /* Obtains one row of nameLists, referring to a specific folder path */
    void nameList(const string& folderPath, vector<string>& nameList);

    /* Erases the contents of a folder, given its path */
    void emptyFolder(const string& folderPath);

    /* Launches a new thread, running the calculations for an element of
    'indexes' */
    void newThread(int whichFolder, int whichFile);


};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void CurveMatching::run() {

    // Saves the current date and time

    auto clock = std::time(nullptr);
    auto tm = *localtime(&clock);

    ostringstream dateStream;
    dateStream << std::put_time(&tm, "%Y-%m-%d");
    date = dateStream.str();

    ostringstream timeStream;
    timeStream << std::put_time(&tm, "%H-%M-%S");
    time = timeStream.str();

    // Removes the files created by previous runs of the program, if any
    emptyFolder("./Graph Data/");

    // Obtains the list of the non-invisible folders and the list of their paths
    // in the Input folder
    listFoldersListFolderPaths("./Input/");

    // Obtains the lists of the names of the files contained in each folder of
    // input data, minus the _exp.txt and _mod.txt appendixes
    nameLists = vector<vector<string>>(folderList.size());
    for (int a=0; a<folderList.size(); ++a)
        nameList(folderPathList[a],nameLists[a]);

    // Generates the Pascal's triangle
    pascalsTriangle = vector<vector<double>>(m,vector<double>(1,1));
    pascalsTriangle[1].push_back(1);
    if (g > 1)
        for (int a=2; a<m; ++a) {
            for (int b=0; b<a-1; ++b)
                pascalsTriangle[a].push_back(
                    pascalsTriangle[a-1][b] + pascalsTriangle[a-1][b+1]);
            pascalsTriangle[a].push_back(1);
        }

    // Initializes 'indexes'
    indexes = vector<vector<Indexes>>(nameLists.size());
    for (int a=0; a<nameLists.size(); ++a)
        indexes[a] = vector<Indexes>(nameLists[a].size());

    // Runs the calculations for each element of indexes
    auto t = vector<vector<thread>>(nameLists.size());
    for (int a=0; a<nameLists.size(); ++a)
        t[a] = vector<thread>(nameLists[a].size());
    for (int a=0; a<nameLists.size(); ++a)
        for (int b=0; b<nameLists[a].size(); ++b)
            t[a][b] = thread(&CurveMatching::newThread,this,a,b);
    for (int a=0; a<nameLists.size(); ++a)
        for (int b=0; b<nameLists[a].size(); ++b)
            t[a][b].join();

}


void CurveMatching::listFoldersListFolderPaths(const string& path) {

    // Obtains the list of the non-invisible folders at 'path'

    const char* pathChar = path.c_str();

    DIR *pdir;
    struct dirent *pent;

    pdir=opendir(pathChar);
    while (pent=readdir(pdir))
        if (pent->d_type == DT_DIR)
            if (pent->d_name[0] != '.')
                folderList.push_back(pent->d_name);
    closedir(pdir);

    // Obtains the list of the paths to the non-invisible folders at 'path'

    string pathString = path;
    if (pathString.back() != '/') pathString += "/";

    for (int a=0; a<folderList.size(); ++a)
        folderPathList.push_back(pathString+folderList[a]+"/");

}



void CurveMatching::nameList(const string& folderPath,
                             vector<string>& nameList) {

    const char* pathChar = folderPath.c_str();

    DIR *pdir;
    struct dirent *pent;

    pdir=opendir(pathChar);
    while (pent=readdir(pdir))
        if (pent->d_type != DT_DIR)
            if (pent->d_name[0] != '.') {
                string fileName = pent->d_name;
                int size = fileName.size();
                if (fileName.back() != '~')
                    if (size > 8)
                        if (fileName.substr(size-8,8) == "_exp.txt")
                            nameList.push_back(fileName.substr(0,size-8));
            }
    closedir(pdir);

}

void CurveMatching::emptyFolder(const string& folderPath) {

    const char* pathChar = folderPath.c_str();

    vector<string> listOfFilesInFolder;

    DIR *pdir;
    struct dirent *pent;

    pdir=opendir(pathChar);
    while (pent=readdir(pdir))
        if (pent->d_type != DT_DIR)
            listOfFilesInFolder.push_back(pent->d_name);
    closedir(pdir);

    for (int a=0; a<listOfFilesInFolder.size(); ++a)
        remove((folderPath+listOfFilesInFolder[a]).c_str());

}

void CurveMatching::newThread(int whichFolder, int whichFile) {

    indexes[whichFolder][whichFile].solve(folderPathList[whichFolder],
                                          folderList[whichFolder],
                                          nameLists[whichFolder][whichFile]);

}



