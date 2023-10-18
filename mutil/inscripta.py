import pandas as pd


def get_designer_input(aledb_uniq_mut_df):
    onyx_example_df = pd.read_csv('./data/DesignTemplateCSV.csv', header=2)
    inscripta_designer_input_var_df = pd.DataFrame(columns=onyx_example_df.columns)

    for _, var in aledb_uniq_mut_df.iterrows():
        assert (len(
            var["clean genes"]) == 1)  # Haven't yet written logic to deal with mutations on more than one features

        rename_d = {
            'SNP': 'Substitution',
            'DEL': 'Deletion',
            'INS': 'Insertion',
        }

        targ_type = ''
        targ_name = ''
        pos_val = -1
        if var['rel nuc pos'] != '':
            pos_val = var['rel nuc pos']
            targ_type = 'Specific CDS'
            targ_name = var['Gene']
        else:
            pos_val = var['abs nuc pos']
            targ_type = 'Chromosome'
            targ_name = 'I'

        inscripta_designer_input_var_df = inscripta_designer_input_var_df.append({
            'INSCRIPTA_EditType': rename_d[var["Mutation Type"]],
            'INSCRIPTA_TargetType': targ_type,
            'INSCRIPTA_TargetName': targ_name,
            'INSCRIPTA_PositionType': 'Specific Position',
            'INSCRIPTA_CoordinateType': 'Nucleotide',  # always Nucleotide for ALEdb mutations
            'INSCRIPTA_PositionValue': pos_val,
            'INSCRIPTA_NumberOfCoordinatesToDelete': var['DEL size'],
            'INSCRIPTA_InsertionSequence': var['SUB or INS seq'],
        }, ignore_index=True)

    # The below tests are duplicated with the pangenome implementation.
    # TODO: find a way to combine (probably as a function in an external file.
    inscripta_designer_input_var_df.fillna('', inplace=True)

    return inscripta_designer_input_var_df
