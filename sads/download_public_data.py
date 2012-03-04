from retriever import VERSION, SCRIPT_LIST, ENGINE_LIST
from retriever.lib.tools import choose_engine, get_opts

public_data = ['BBS', 'FIA', 'MCDB']

for dataset in public_data:
    script_list = SCRIPT_LIST()
    opts = get_opts(script_list, args=['install', dataset, '-e', 's', '-f',
                                       'downloaded_data.sqlite'])
    script = opts["script"]
    engine = choose_engine(opts)
    if isinstance(script, list):
        for dataset in script:
            print "=> Installing", dataset.name
            dataset.download(engine, debug=debug)
    else:
        script.download(engine)
print "Datasets successfully downloaded."