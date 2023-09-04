# lib/tasks/csv_load.rake

namespace :csv_load do
    desc 'Load compound properties from PubChem API and save as CSV files'
    task :load_properties => :environment do
      # Replace the following lines with your logic to fetch the compounds from your database or any other source
      
     inhibitors = Inhibitor.all
  
      inhibitors.each do |inhibitor|
        puts "Loading properties for #{inhibitor[:name]}..."
        PubChemService.get_compound_properties(inhibitor[:smiles], inhibitor[:name])
        puts "Properties for #{inhibitor[:name]} loaded and saved as CSV."
      end
    end
  end
  