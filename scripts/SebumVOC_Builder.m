%% ========================================================================
%  SEBUM VOC DATABASE BUILDER v3.1 (fixed)
%  Input:  XLSX file (compounds in column A)
%  Output: Same or new XLSX file with all database columns added
%  APIs:   PubChem PUG REST, HMDB, KEGG, ChEBI, ClassyFire, NIST
%% ========================================================================
function SebumVOC_Builder()
clc;
fprintf('==============================================\n');
fprintf('   SEBUM VOC DATABASE BUILDER v3.1\n');
fprintf('==============================================\n\n');

%% ---- STEP 1: Select input XLSX file ------------------------------------
[fname, fpath] = uigetfile('*.xlsx', 'Select input XLSX file with compound names');
if isequal(fname, 0)
    fprintf('No file selected. Exiting.\n');
    return;
end
inputFile = fullfile(fpath, fname);
fprintf('Input file: %s\n', inputFile);

%% ---- STEP 2: Read compounds from column A ------------------------------
try
    rawData = readcell(inputFile, 'Sheet', 1);
catch ME
    error('Cannot read Excel file: %s', ME.message);
end

compounds = {};
for r = 1:size(rawData, 1)
    val = rawData{r, 1};
    if ischar(val) || isstring(val)
        val = strtrim(char(val));
        if ~isempty(val)
            compounds{end+1} = val; %#ok<AGROW>
        end
    elseif isnumeric(val) && ~isnan(val)
        compounds{end+1} = num2str(val); %#ok<AGROW>
    end
end
nCompounds = numel(compounds);
if nCompounds == 0
    error('No compounds found in column A of the selected file.');
end
fprintf('Found %d compounds.\n\n', nCompounds);

%% ---- STEP 3: Configuration ---------------------------------------------
cfg.delay      = 0.5;
cfg.maxRetry   = 3;
cfg.verbose    = true;
cfg.saveEvery  = 20;

%% ---- STEP 4: Define output columns -------------------------------------
headers = {
    'Input_Name'
    'Preferred_IUPAC_Name'
    'Common_Synonyms'
    'CAS_Number'
    'PubChem_CID'
    'InChIKey'
    'InChI'
    'Canonical_SMILES'
    'Molecular_Formula'
    'Molecular_Weight_g_mol'
    'RI_Polar_Column_NIST'
    'RI_Nonpolar_Column_NIST'
    'NIST_WebBook_URL'
    'ChEBI_ID'
    'ChEBI_URL'
    'HMDB_ID'
    'HMDB_URL'
    'KEGG_Compound_ID'
    'KEGG_URL'
    'PubChem_URL'
    'ClassyFire_Kingdom'
    'ClassyFire_Superclass'
    'ClassyFire_Class'
    'ClassyFire_Subclass'
    'Metabolic_Pathways'
    'KEGG_Pathway_IDs'
    'Disease_Associations'
    'Synthesized_in_Disease'
    'Biosynthesis_Info'
    'Notes'
    };
nCols = numel(headers);

%% ---- STEP 5: Initialise results matrix ---------------------------------
results = repmat({'Not found'}, nCompounds, nCols);

%% ---- STEP 6: Main processing loop --------------------------------------
fprintf('Starting database queries...\n');
fprintf('----------------------------------------------\n');

for i = 1:nCompounds
    rawName = strtrim(compounds{i});
    if isempty(rawName), continue; end
    fprintf('[%d/%d] %s\n', i, nCompounds, rawName);

    results{i, getCol('Input_Name', headers)} = rawName;
    qName = cleanName(rawName);

    %% --- PubChem --------------------------------------------------------
    [cid, pc] = pubchemQuery(qName, cfg);
    if ~isempty(cid)
        results{i, getCol('PubChem_CID',            headers)} = cid;
        results{i, getCol('Preferred_IUPAC_Name',   headers)} = pc.iupac;
        results{i, getCol('InChIKey',               headers)} = pc.inchikey;
        results{i, getCol('InChI',                  headers)} = pc.inchi;
        results{i, getCol('Canonical_SMILES',       headers)} = pc.smiles;
        results{i, getCol('Molecular_Formula',      headers)} = pc.formula;
        results{i, getCol('Molecular_Weight_g_mol', headers)} = pc.mw;
        results{i, getCol('PubChem_URL',            headers)} = ...
            ['https://pubchem.ncbi.nlm.nih.gov/compound/' cid];

        [syns, cas] = pubchemSynonyms(cid, cfg);
        results{i, getCol('Common_Synonyms', headers)} = syns;
        results{i, getCol('CAS_Number',      headers)} = cas;

        if ~strcmp(pc.inchikey, 'Not found')
            cf = classyfireQuery(pc.inchikey, cfg);
            results{i, getCol('ClassyFire_Kingdom',    headers)} = cf.kingdom;
            results{i, getCol('ClassyFire_Superclass', headers)} = cf.superclass;
            results{i, getCol('ClassyFire_Class',      headers)} = cf.class;
            results{i, getCol('ClassyFire_Subclass',   headers)} = cf.subclass;
        end
    else
        if cfg.verbose
            fprintf('  [!] PubChem: not found\n');
        end
        results{i, getCol('Notes', headers)} = 'Not found in PubChem - verify name';
    end

    %% --- NIST WebBook RI ------------------------------------------------
    nistURL = ['https://webbook.nist.gov/cgi/cbook.cgi?Name=' ...
               urlEncode(qName) '&Units=SI&cGC=on'];
    results{i, getCol('NIST_WebBook_URL', headers)} = nistURL;
    [riP, riNP] = nistRIquery(qName, cfg);
    results{i, getCol('RI_Polar_Column_NIST',    headers)} = riP;
    results{i, getCol('RI_Nonpolar_Column_NIST', headers)} = riNP;

    %% --- ChEBI ----------------------------------------------------------
    chebi = chebiQuery(qName, cid, cfg);
    if ~isempty(chebi.id)
        results{i, getCol('ChEBI_ID',  headers)} = chebi.id;
        results{i, getCol('ChEBI_URL', headers)} = ...
            ['https://www.ebi.ac.uk/chebi/searchId.do?chebiId=' chebi.id];
    end

    %% --- HMDB -----------------------------------------------------------
    hmdb = hmdbQuery(qName, cid, cfg);
    if ~isempty(hmdb.id)
        results{i, getCol('HMDB_ID',               headers)} = hmdb.id;
        results{i, getCol('HMDB_URL',              headers)} = ...
            ['https://hmdb.ca/metabolites/' hmdb.id];
        results{i, getCol('Metabolic_Pathways',    headers)} = hmdb.pathways;
        results{i, getCol('Disease_Associations',  headers)} = hmdb.diseases;
        results{i, getCol('Synthesized_in_Disease',headers)} = hmdb.synthesized;
        results{i, getCol('Biosynthesis_Info',     headers)} = hmdb.biosynth;
    end

    %% --- KEGG -----------------------------------------------------------
    kegg = keggQuery(qName, cfg);
    if ~isempty(kegg.id)
        results{i, getCol('KEGG_Compound_ID', headers)} = kegg.id;
        results{i, getCol('KEGG_URL',         headers)} = ...
            ['https://www.genome.jp/entry/' kegg.id];
        results{i, getCol('KEGG_Pathway_IDs', headers)} = kegg.pathways;
        if strcmp(results{i, getCol('Metabolic_Pathways', headers)}, 'Not found')
            results{i, getCol('Metabolic_Pathways', headers)} = kegg.pathways;
        end
    end

    %% --- Quality flag ---------------------------------------------------
    if strcmp(results{i, getCol('Notes', headers)}, 'Not found')
        results{i, getCol('Notes', headers)} = 'OK';
    end

    pause(cfg.delay);

    %% --- Auto-save ------------------------------------------------------
    if mod(i, cfg.saveEvery) == 0
        tmpFile = fullfile(fpath, ['TEMP_' fname]);
        writeExcel(tmpFile, headers, results, i);
        fprintf('  --> Auto-saved %d rows to %s\n', i, tmpFile);
    end
end

%% ---- STEP 7: Write final output ----------------------------------------
[~, baseName, ~] = fileparts(fname);
outputFile = fullfile(fpath, [baseName '_DATABASE.xlsx']);
writeExcel(outputFile, headers, results, nCompounds);
fprintf('\n==============================================\n');
fprintf('COMPLETE! Saved to:\n%s\n', outputFile);
fprintf('==============================================\n');
end

%% ========================================================================
%%  HELPER: Column index lookup
%% ========================================================================
function idx = getCol(name, headers)
    idx = find(strcmp(headers, name), 1);
    if isempty(idx)
        error('Column not found: %s', name);
    end
end

%% ========================================================================
%%  HELPER: Clean compound name for API queries
%% ========================================================================
function out = cleanName(name)
    out = name;
    out = strrep(out, '.alpha.', 'alpha');
    out = strrep(out, '.beta.',  'beta');
    out = strrep(out, '.gamma.', 'gamma');
    out = strrep(out, '.delta.', 'delta');
    out = strrep(out, char(945), 'alpha');
    out = strrep(out, char(946), 'beta');
    out = strrep(out, char(947), 'gamma');
    out = strrep(out, '(.+/-.)', '');
    out = strrep(out, '[1R-',   '');
    out = strrep(out, '[3R-',   '');
    out = strrep(out, '(1R)',   '');
    out = strrep(out, '(E)-',   '');
    out = strrep(out, '(Z)-',   '');
    out = regexprep(out, '\s+', ' ');
    out = strtrim(out);
end

%% ========================================================================
%%  HELPER: URL encode
%% ========================================================================
function encoded = urlEncode(str)
    try
        encoded = urlencode(str);
    catch
        encoded = str;
        encoded = strrep(encoded, ' ',  '+');
        encoded = strrep(encoded, ',',  '%2C');
        encoded = strrep(encoded, '(',  '%28');
        encoded = strrep(encoded, ')',  '%29');
        encoded = strrep(encoded, '[',  '%5B');
        encoded = strrep(encoded, ']',  '%5D');
        encoded = strrep(encoded, '''', '%27');
        encoded = strrep(encoded, '#',  '%23');
    end
end

%% ========================================================================
%%  OUTPUT: Write results to Excel
%% ========================================================================
function writeExcel(outputFile, headers, results, nRows)
    outCell = [headers'; results(1:nRows, :)];
    try
        writecell(outCell, outputFile, 'Sheet', 1);
    catch ME
        warning('Could not write to %s: %s', outputFile, ME.message);
        try
            xlswrite(outputFile, outCell, 1); %#ok<XLSWT>
        catch ME2
            warning('xlswrite also failed: %s', ME2.message);
        end
    end
end

%% ========================================================================
%%  API: PubChem – main compound query (with fallback)
%% ========================================================================
function [cid, pc] = pubchemQuery(compoundName, cfg)
    cid = '';
    pc  = struct('iupac','Not found','inchikey','Not found','inchi','Not found',...
                 'smiles','Not found','formula','Not found','mw','Not found');

    cidURL = ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' ...
              urlEncode(compoundName) '/cids/JSON'];
    cidJSON = robustWebRead(cidURL, cfg);

    % Fallback: PubChem autocomplete
    if isempty(cidJSON)
        acURL = ['https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/' ...
                 urlEncode(compoundName) '/JSON?limit=1'];
        acJSON = robustWebRead(acURL, cfg);
        if ~isempty(acJSON)
            try
                acData = jsondecode(acJSON);
                if isfield(acData,'dictionary_terms') && ...
                   isfield(acData.dictionary_terms,'compound')
                    altName = acData.dictionary_terms.compound;
                    if iscell(altName), altName = altName{1}; end
                    altName = char(altName);
                    cidURL2 = ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' ...
                               urlEncode(altName) '/cids/JSON'];
                    cidJSON = robustWebRead(cidURL2, cfg);
                end
            catch
            end
        end
    end

    if isempty(cidJSON), return; end

    try
        cidData = jsondecode(cidJSON);
    catch
        return;
    end

    if isfield(cidData,'IdentifierList') && isfield(cidData.IdentifierList,'CID')
        cidList = cidData.IdentifierList.CID;
        if isnumeric(cidList)
            cid = num2str(cidList(1));
        end
    end
    if isempty(cid), return; end

    propList = 'IUPACName,InChIKey,InChI,CanonicalSMILES,MolecularFormula,MolecularWeight';
    propURL  = ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' ...
                cid '/property/' propList '/JSON'];
    propJSON = robustWebRead(propURL, cfg);
    if isempty(propJSON), return; end

    try
        data = jsondecode(propJSON);
    catch
        return;
    end

    if isfield(data,'PropertyTable') && isfield(data.PropertyTable,'Properties')
        props = data.PropertyTable.Properties;
        if iscell(props), props = props{1}; end
        if isfield(props,'IUPACName'),        pc.iupac    = char(props.IUPACName);       end
        if isfield(props,'InChIKey'),         pc.inchikey = char(props.InChIKey);         end
        if isfield(props,'InChI'),            pc.inchi    = char(props.InChI);            end
        if isfield(props,'CanonicalSMILES'),  pc.smiles   = char(props.CanonicalSMILES); end
        if isfield(props,'MolecularFormula'), pc.formula  = char(props.MolecularFormula); end
        if isfield(props,'MolecularWeight')
            if isnumeric(props.MolecularWeight)
                pc.mw = num2str(props.MolecularWeight);
            else
                pc.mw = char(props.MolecularWeight);
            end
        end
    end
end

%% ========================================================================
%%  API: PubChem – Synonyms & CAS extraction
%% ========================================================================
function [synStr, casStr] = pubchemSynonyms(cid, cfg)
    synStr = 'Not found';
    casStr = 'Not found';

    url = ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' ...
           cid '/synonyms/JSON'];
    jsonText = robustWebRead(url, cfg);
    if isempty(jsonText), return; end

    try
        data = jsondecode(jsonText);
    catch
        return;
    end

    if isfield(data,'InformationList') && isfield(data.InformationList,'Information')
        info = data.InformationList.Information;
        if iscell(info), info = info{1}; end
        if isfield(info,'Synonym')
            syns = info.Synonym;
            if iscell(syns)
                syns = cellfun(@char, syns, 'UniformOutput', false);
            elseif isstring(syns)
                syns = cellstr(syns);
            else
                syns = {char(syns)};
            end
            casPattern = '^\d{2,7}-\d{2}-\d$';
            for s = 1:numel(syns)
                if ~isempty(regexp(syns{s}, casPattern, 'once'))
                    casStr = syns{s};
                    break;
                end
            end
            maxSyn = min(10, numel(syns));
            synStr = strjoin(syns(1:maxSyn), '; ');
        end
    end
end

%% ========================================================================
%%  API: ClassyFire – Taxonomy via InChIKey
%% ========================================================================
function cf = classyfireQuery(inchikey, cfg)
    cf = struct('kingdom','Not found','superclass','Not found',...
                'class','Not found','subclass','Not found');

    url = ['http://classyfire.wishartlab.com/entities/' char(inchikey) '.json'];
    jsonText = robustWebRead(url, cfg);

    if isempty(jsonText)
        url2 = ['https://gnps-classyfire.ucsd.edu/entities/' char(inchikey) '.json'];
        jsonText = robustWebRead(url2, cfg);
    end
    if isempty(jsonText), return; end

    try
        data = jsondecode(jsonText);
    catch
        return;
    end

    if isfield(data,'kingdom')    && isstruct(data.kingdom)    && isfield(data.kingdom,'name'),    cf.kingdom    = char(data.kingdom.name);    end
    if isfield(data,'superclass') && isstruct(data.superclass) && isfield(data.superclass,'name'), cf.superclass = char(data.superclass.name); end
    if isfield(data,'class')      && isstruct(data.class)      && isfield(data.class,'name'),      cf.class      = char(data.class.name);      end
    if isfield(data,'subclass')   && isstruct(data.subclass)   && isfield(data.subclass,'name'),   cf.subclass   = char(data.subclass.name);   end
end

%% ========================================================================
%%  API: NIST WebBook – Retention Index scraping
%% ========================================================================
function [riPolar, riNonpolar] = nistRIquery(compoundName, cfg)
    riPolar    = 'Not found';
    riNonpolar = 'Not found';

    url = ['https://webbook.nist.gov/cgi/cbook.cgi?Name=' ...
           urlEncode(compoundName) '&Units=SI&cGC=on'];
    htmlText = robustWebRead(url, cfg);
    if isempty(htmlText), return; end

    npPatterns = {'Standard\s+non-polar.*?(\d{3,5})', ...
                  'Non-polar.*?(\d{3,5})', ...
                  'SE-30.*?(\d{3,5})', 'OV-1.*?(\d{3,5})', ...
                  'DB-1.*?(\d{3,5})', 'HP-5.*?(\d{3,5})', 'DB-5.*?(\d{3,5})'};
    for p = 1:numel(npPatterns)
        tok = regexp(htmlText, npPatterns{p}, 'tokens', 'once');
        if ~isempty(tok), riNonpolar = tok{1}; break; end
    end

    pPatterns = {'Standard\s+polar.*?(\d{3,5})', ...
                 'Polar.*?(\d{3,5})', ...
                 'Carbowax.*?(\d{3,5})', 'DB-WAX.*?(\d{3,5})', ...
                 'FFAP.*?(\d{3,5})', 'PEG.*?(\d{3,5})'};
    for p = 1:numel(pPatterns)
        tok = regexp(htmlText, pPatterns{p}, 'tokens', 'once');
        if ~isempty(tok), riPolar = tok{1}; break; end
    end

    if cfg.verbose
        fprintf('  [NIST] RI Polar=%s  Non-polar=%s\n', riPolar, riNonpolar);
    end
end

%% ========================================================================
%%  API: ChEBI – via PubChem cross-refs + OLS4 REST (3 arguments!)
%% ========================================================================
function chebi = chebiQuery(compoundName, pubchemCID, cfg)
    chebi = struct('id', '');

    % Strategy 1: PubChem cross-references (most reliable)
    if ~isempty(pubchemCID)
        xrefURL = ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' ...
                   pubchemCID '/xrefs/RegistryID/JSON'];
        xrefJSON = robustWebRead(xrefURL, cfg);
        if ~isempty(xrefJSON)
            try
                xrefData = jsondecode(xrefJSON);
                if isfield(xrefData,'InformationList') && ...
                   isfield(xrefData.InformationList,'Information')
                    info = xrefData.InformationList.Information;
                    if iscell(info), info = info{1}; end
                    if isfield(info,'RegistryID')
                        regIDs = info.RegistryID;
                        if iscell(regIDs)
                            regIDs = cellfun(@char, regIDs, 'UniformOutput', false);
                        elseif isstring(regIDs)
                            regIDs = cellstr(regIDs);
                        else
                            regIDs = {char(regIDs)};
                        end
                        for r = 1:numel(regIDs)
                            if contains(regIDs{r}, 'CHEBI')
                                chebi.id = regIDs{r};
                                if ~contains(chebi.id, ':')
                                    chebi.id = strrep(chebi.id, 'CHEBI', 'CHEBI:');
                                end
                                if cfg.verbose
                                    fprintf('  [ChEBI] %s (via PubChem xref)\n', chebi.id);
                                end
                                return;
                            end
                        end
                    end
                end
            catch
            end
        end
    end

    % Strategy 2: OLS4 REST API exact search
    url1 = ['https://www.ebi.ac.uk/ols4/api/search?q=' ...
            urlEncode(compoundName) '&ontology=chebi&exact=true&rows=1'];
    jsonText = robustWebRead(url1, cfg);
    if ~isempty(jsonText)
        chebi.id = parseOLSforChEBI(jsonText);
        if ~isempty(chebi.id)
            if cfg.verbose, fprintf('  [ChEBI] %s (via OLS4 exact)\n', chebi.id); end
            return;
        end
    end

    % Strategy 3: OLS4 non-exact
    url2 = ['https://www.ebi.ac.uk/ols4/api/search?q=' ...
            urlEncode(compoundName) '&ontology=chebi&rows=1'];
    jsonText = robustWebRead(url2, cfg);
    if ~isempty(jsonText)
        chebi.id = parseOLSforChEBI(jsonText);
        if ~isempty(chebi.id)
            if cfg.verbose, fprintf('  [ChEBI] %s (via OLS4)\n', chebi.id); end
            return;
        end
    end

    if cfg.verbose
        fprintf('  [ChEBI] Not found\n');
    end
end

function chebiID = parseOLSforChEBI(jsonText)
    chebiID = '';
    try
        data = jsondecode(jsonText);
    catch
        return;
    end
    if isfield(data,'response') && isfield(data.response,'docs')
        docs = data.response.docs;
        if iscell(docs) && ~isempty(docs)
            doc = docs{1};
        elseif isstruct(docs) && ~isempty(docs)
            doc = docs(1);
        else
            return;
        end
        if isfield(doc,'obo_id')
            chebiID = char(doc.obo_id);
        elseif isfield(doc,'short_form')
            chebiID = strrep(char(doc.short_form), '_', ':');
        end
    end
end

%% ========================================================================
%%  API: HMDB – via PubChem cross-refs + HMDB XML (3 arguments!)
%% ========================================================================
function hmdb = hmdbQuery(compoundName, pubchemCID, cfg)
    hmdb = struct('id','','pathways','Not found','diseases','Not found',...
                  'synthesized','Not found','biosynth','Not found');

    hmdbID = '';

    % Strategy 1: PubChem cross-references
    if ~isempty(pubchemCID)
        xrefURL = ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' ...
                   pubchemCID '/xrefs/RegistryID/JSON'];
        xrefJSON = robustWebRead(xrefURL, cfg);
        if ~isempty(xrefJSON)
            try
                xrefData = jsondecode(xrefJSON);
                if isfield(xrefData,'InformationList') && ...
                   isfield(xrefData.InformationList,'Information')
                    info = xrefData.InformationList.Information;
                    if iscell(info), info = info{1}; end
                    if isfield(info,'RegistryID')
                        regIDs = info.RegistryID;
                        if iscell(regIDs)
                            regIDs = cellfun(@char, regIDs, 'UniformOutput', false);
                        elseif isstring(regIDs)
                            regIDs = cellstr(regIDs);
                        else
                            regIDs = {char(regIDs)};
                        end
                        for r = 1:numel(regIDs)
                            if startsWith(regIDs{r}, 'HMDB')
                                hmdbID = regIDs{r};
                                break;
                            end
                        end
                    end
                end
            catch
            end
        end
    end

    % Strategy 2: PubChem synonyms
    if isempty(hmdbID) && ~isempty(pubchemCID)
        synURL = ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' ...
                  pubchemCID '/synonyms/JSON'];
        synJSON = robustWebRead(synURL, cfg);
        if ~isempty(synJSON)
            try
                synData = jsondecode(synJSON);
                if isfield(synData,'InformationList') && ...
                   isfield(synData.InformationList,'Information')
                    sInfo = synData.InformationList.Information;
                    if iscell(sInfo), sInfo = sInfo{1}; end
                    if isfield(sInfo,'Synonym')
                        allSyns = sInfo.Synonym;
                        if iscell(allSyns)
                            allSyns = cellfun(@char, allSyns, 'UniformOutput', false);
                        else
                            allSyns = {char(allSyns)};
                        end
                        for ss = 1:numel(allSyns)
                            tok = regexp(allSyns{ss}, '(HMDB\d{4,7})', 'tokens', 'once');
                            if ~isempty(tok)
                                hmdbID = tok{1};
                                break;
                            end
                        end
                    end
                end
            catch
            end
        end
    end

    if isempty(hmdbID)
        if cfg.verbose, fprintf('  [HMDB] Not found\n'); end
        return;
    end

    % Normalize HMDB ID to 7-digit format
    hmdbNum = regexp(hmdbID, 'HMDB(\d+)', 'tokens', 'once');
    if ~isempty(hmdbNum)
        numStr = hmdbNum{1};
        if length(numStr) < 7
            numStr = [repmat('0', 1, 7-length(numStr)) numStr];
        end
        hmdbID = ['HMDB' numStr];
    end

    hmdb.id = hmdbID;
    if cfg.verbose, fprintf('  [HMDB] %s\n', hmdb.id); end

    % Fetch detailed XML from HMDB
    xmlURL = ['https://hmdb.ca/metabolites/' hmdbID '.xml'];
    xmlText = robustWebRead(xmlURL, cfg);
    if isempty(xmlText), return; end

    % Extract pathways
    pwPattern = '<pathway>.*?<name>(.*?)</name>.*?</pathway>';
    pwTokens = regexp(xmlText, pwPattern, 'tokens');
    if ~isempty(pwTokens)
        pwNames = cellfun(@(c) c{1}, pwTokens, 'UniformOutput', false);
        hmdb.pathways = strjoin(unique(pwNames), '; ');
    end

    % Extract diseases
    disPattern = '<disease>.*?<name>(.*?)</name>.*?</disease>';
    disTokens = regexp(xmlText, disPattern, 'tokens');
    if ~isempty(disTokens)
        disNames = cellfun(@(c) c{1}, disTokens, 'UniformOutput', false);
        hmdb.diseases = strjoin(unique(disNames), '; ');
    end

    % Extract description for disease-related synthesis
    descPattern = '<description>(.*?)</description>';
    descTokens = regexp(xmlText, descPattern, 'tokens', 'once');
    if ~isempty(descTokens)
        description = descTokens{1};
        description = regexprep(description, '<[^>]+>', '');
        description = strtrim(description);
        diseaseKW = {'cancer','diabetes','disease','disorder','syndrome',...
                     'infection','inflammation','tumor','tumour','obesity',...
                     'parkinson','alzheimer','metabolic','biomarker'};
        for dk = 1:numel(diseaseKW)
            if ~isempty(regexpi(description, diseaseKW{dk}))
                if length(description) > 500
                    description = [description(1:497) '...']; %#ok<AGROW>
                end
                hmdb.synthesized = description;
                break;
            end
        end
    end

    % Extract biosynthesis references
    bioPattern = '<synthesis_reference>(.*?)</synthesis_reference>';
    bioTokens = regexp(xmlText, bioPattern, 'tokens');
    if ~isempty(bioTokens)
        bioInfo = cellfun(@(c) strtrim(c{1}), bioTokens, 'UniformOutput', false);
        hmdb.biosynth = strjoin(bioInfo, '; ');
    end
end

%% ========================================================================
%%  API: KEGG – Search by compound name (fixed URL format)
%% ========================================================================
function kegg = keggQuery(compoundName, cfg)
    kegg = struct('id','','pathways','Not found');

    findURL = ['https://rest.kegg.jp/find/compound/' urlEncode(compoundName)];
    findText = robustWebRead(findURL, cfg);
    if isempty(findText), return; end

    lines = strsplit(findText, '\n');
    if isempty(lines) || isempty(strtrim(lines{1})), return; end

    firstLine = strtrim(lines{1});
    parts = strsplit(firstLine, '\t');
    if isempty(parts), return; end

    keggEntry = strtrim(parts{1});
    keggEntry = strrep(keggEntry, 'cpd:', '');
    kegg.id = keggEntry;

    if cfg.verbose, fprintf('  [KEGG] %s\n', kegg.id); end

    % FIXED URL: https://rest.kegg.jp/link/pathway/C00084
    % (NOT /link/pathway/compound/C00084)
    linkURL = ['https://rest.kegg.jp/link/pathway/' kegg.id];
    linkText = robustWebRead(linkURL, cfg);
    if isempty(linkText), return; end

    pwLines = strsplit(linkText, '\n');
    pwIDs = {};
    for p = 1:numel(pwLines)
        ln = strtrim(pwLines{p});
        if isempty(ln), continue; end
        cols = strsplit(ln, '\t');
        if numel(cols) >= 2
            pwID = strrep(strtrim(cols{2}), 'path:', '');
            pwID = regexprep(pwID, '^[a-z]{2,4}(?=\d)', 'map');
            pwIDs{end+1} = pwID; %#ok<AGROW>
        end
    end

    if isempty(pwIDs), return; end

    % Get pathway names (limit to first 8)
    uniquePW = unique(pwIDs);
    maxPW = min(8, numel(uniquePW));
    pwNamesAll = {};
    for pp = 1:maxPW
        pwGetURL = ['https://rest.kegg.jp/get/' uniquePW{pp}];
        pwGetText = robustWebRead(pwGetURL, cfg);
        if ~isempty(pwGetText)
            nameMatch = regexp(pwGetText, 'NAME\s+(.*?)(?:\n|$)', 'tokens', 'once');
            if ~isempty(nameMatch)
                pwName = strtrim(nameMatch{1});
                pwName = regexprep(pwName, '\s*-\s*Reference pathway.*', '');
                pwName = regexprep(pwName, '\s*-\s*Homo sapiens.*', '');
                pwNamesAll{end+1} = [uniquePW{pp} ': ' pwName]; %#ok<AGROW>
            else
                pwNamesAll{end+1} = uniquePW{pp}; %#ok<AGROW>
            end
        end
        pause(cfg.delay * 0.5);
    end

    if ~isempty(pwNamesAll)
        kegg.pathways = strjoin(pwNamesAll, '; ');
    else
        kegg.pathways = strjoin(uniquePW, '; ');
    end
end

%% ========================================================================
%%  HELPER: Robust web read with retries
%% ========================================================================
function responseText = robustWebRead(url, cfg)
    responseText = '';
    for attempt = 1:cfg.maxRetry
        try
            options = weboptions('Timeout', 30, ...
                                 'ContentType', 'text', ...
                                 'CharacterEncoding', 'UTF-8');
            options.HeaderFields = {'User-Agent', ...
                'MATLAB SebumVOC_Builder/3.1 (research use)'};
            responseText = webread(url, options);
            return;
        catch ME
            errMsg = ME.message;
            is404 = contains(errMsg,'404') || contains(errMsg,'NotFound');
            is400 = contains(errMsg,'400') || contains(errMsg,'Bad Request');
            is403 = contains(errMsg,'403') || contains(errMsg,'Forbidden');
            is500 = contains(errMsg,'500') || contains(errMsg,'Internal Server');

            if is404 || is400
                if cfg.verbose
                    fprintf('  [HTTP %s] %s\n', extractHTTPcode(errMsg), shortenURL(url));
                end
                return;
            end
            if (is403 || is500) && attempt >= 2
                if cfg.verbose
                    fprintf('  [HTTP %s] %s\n', extractHTTPcode(errMsg), shortenURL(url));
                end
                return;
            end
            if cfg.verbose && attempt == 1
                fprintf('  [WARN] attempt %d/%d: %s\n', attempt, cfg.maxRetry, shortenURL(url));
            end
            if attempt < cfg.maxRetry
                pause(cfg.delay * attempt);
            end
        end
    end
    if cfg.verbose
        fprintf('  [ERR] All %d retries failed: %s\n', cfg.maxRetry, shortenURL(url));
    end
end

%% ========================================================================
%%  HELPER: Extract HTTP status code from error message
%% ========================================================================
function code = extractHTTPcode(msg)
    tok = regexp(msg, '(\d{3})', 'tokens', 'once');
    if ~isempty(tok)
        code = tok{1};
    else
        code = '???';
    end
end

%% ========================================================================
%%  HELPER: Shorten URL for display
%% ========================================================================
function s = shortenURL(url)
    if length(url) > 80
        s = [url(1:77) '...'];
    else
        s = url;
    end
end

