pipeline {
    agent { node { label 'linux' } }
    stages {
        stage('store') {
            steps {
                sh 'set -e; if [ ! -r testing/.stored ]; then cd testing; ./store-ref-data.sh; touch .stored; fi'
            }
        }
        stage('test') {
            steps {
                sh 'set -e; cd testing; TERM=dumb ./run-tests.sh'
            }
        }
    }
}
