def get_operon_name_set(RegulonDB_ID, TU_objects_df, TU_df, operon_df):
    op_set = set()
    TU_IDs = TU_objects_df[TU_objects_df["TU_OBJECT_ID"]==RegulonDB_object_ID]["TRANSCRIPTION_UNIT_ID"]
    for TU_ID in TU_IDs:
        op_IDs = TU_df[TU_df["TRANSCRIPTION_UNIT_ID"]==TU_ID]["OPERON_ID"]
        for op_ID in op_IDs:
            op_name = operon_df[operon_df["OPERON_ID"]==op_ID].iloc[0]["OPERON_NAME"]  # Assuming that each operon only has one name
            op_set.add(op_name)
    return op_set


def get_operon_ID_set(RegulonDB_object_ID, TU_objects_df, TU_df, operon_df):
    op_set = set()
    TU_IDs = TU_objects_df[TU_objects_df["TU_OBJECT_ID"]==RegulonDB_object_ID]["TRANSCRIPTION_UNIT_ID"]
    for TU_ID in TU_IDs:
        op_IDs = TU_df[TU_df["TRANSCRIPTION_UNIT_ID"]==TU_ID]["OPERON_ID"]
        for op_ID in op_IDs:
            op_set.add(op_ID)
    return op_set