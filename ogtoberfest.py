import sys
import os

ogtoberfest_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'src'))
sys.path.insert(0, ogtoberfest_dir)

if __name__ == "__main__":
    from ogtoberfest.main import main
    args = sys.argv[1:]
    main(args)
